;; min cx
;; s.t. Ax = b
;;  l < x < u
;;
;; y: dual mults for Ax = b,
;; z, w: dual mults for l < x and x < u

(defstruct (pdas-state
            (:conc-name #:pdas-))
  nvars ncons
  c
  triplets
  %A b
  l u
  x y w z)

(defun pdas-A (state)
  (or (pdas-%A state)
      (let ((A (make-sparse-from-triplet-vector
                (pdas-ncons state)
                (pdas-nvars state)
                (pdas-triplets state))))
        (setf (pdas-%A state) A))))

(defun free-sparse (sparse)
  (with-alien ((sparse (* cholmod-sparse)
                       :local sparse))
    (assert (/= (cholmod-free-sparse (addr sparse)
                                     *cholmod-common*)
                0))))

(defun free-pdas-A (state)
  (let ((sparse (shiftf (pdas-%A state) nil)))
    (when sparse
      (free-sparse sparse))))

(defvar *clamp* 1d8)

(defun clamp (vector)
  (let ((min (- *clamp*))
        (max *clamp*))
    (map '(simple-array double-float 1)
         (lambda (x)
           (max min (min x max)))
         vector)))

(define-modify-macro maxf (y) max)
(define-modify-macro mulf (y) *)

(defun scale-constraints (triplets rhs)
  #+nil(return-from scale-constraints (values triplets rhs))
  (let ((row-sums (make-hash-table))
        (rhs      (copy-seq rhs)))
    (map nil (lambda (triplet)
               (maxf (gethash (triplet-row triplet) row-sums 0d0)
                     (abs (triplet-value triplet))))
         triplets)
    (maphash (lambda (row sum)
               (let ((scale (setf (gethash row row-sums)
                                  (if (< sum 1d-6)
                                      1d0
                                      (/ sum)))))
                 (mulf (aref rhs row) scale)))
             row-sums)
    (values (map 'simple-vector
                 (lambda (triplet)
                   (let ((scale (gethash (triplet-row triplet)
                                         row-sums)))
                     (make-triplet :row (triplet-row triplet)
                                   :col (triplet-col triplet)
                                   :value (* (triplet-value triplet) scale))))
                 triplets)
            rhs)))

(defun make-pdas (sf)
  (let* ((nvars (sf-nvars sf))
         (ncons (sf-ncons sf))
         (c (let ((v (matlisp:make-real-matrix-dim nvars 1)))
              (loop for (i . x) across (sf-c sf)
                    do (setf (matlisp:matrix-ref v i) x))
              v))
         (x (matlisp:make-real-matrix-dim nvars 1))
         (l (matlisp:make-real-matrix
             (clamp (sf-l sf))))
         (u (matlisp:make-real-matrix
             (clamp (sf-u sf))))
         (z (matlisp:make-real-matrix-dim nvars 1))
         (w (matlisp:make-real-matrix-dim nvars 1)))
    (dotimes (i (sf-nvars sf))
      (when (< (- (matlisp:matrix-ref u i)
                  (matlisp:matrix-ref l i))
               1d-6)
        (decf (matlisp:matrix-ref l i) 5d-7)
        (incf (matlisp:matrix-ref u i) 5d7))
      (let* ((l (aref (sf-l sf) i))
             (u (aref (sf-u sf) i))
             (delta (- u l))
             (xi (cond ((and (< l -1d10) (> u 1d10))
                        0d0)
                       ((< l -1d6)
                        (- u (min (/ delta 2)
                                  (1+ (* (abs u) .1d0)))))
                       ((> u 1d6)
                        (+ l (min (/ delta 2)
                                  (1+ (* (abs l) .1d0)))))
                       (t
                        (/ (+ l u) 2)))))
        (setf (matlisp:matrix-ref x i) xi)))
    (dotimes (i nvars)
      (let ((c (matlisp:matrix-ref c i)))
        (setf (values (matlisp:matrix-ref z i)
                      (matlisp:matrix-ref w i))
              (cond ((or (zerop c))
                     (values 1d0 1d0))
                    ((minusp c)
                     (values 1d0 (+ 1d0 (- c))))
                    (t
                     (values (+ 1d0 c) 1d0))))))
    (multiple-value-bind (triplets rhs)
        (scale-constraints (sf-A sf) (sf-b sf))
      (make-pdas-state
       :nvars nvars
       :ncons ncons
       :c c
       :triplets triplets
       :%A nil
       :b (matlisp:make-real-matrix rhs)
       :l l
       :u u
       :x x
       :y (matlisp:zeros ncons 1)
       :z z
       :w w))))

(defun violation (state)
  (let* ((x (pdas-x state))
         (l (matlisp:m- x (pdas-l state)))
         (u (matlisp:m- (pdas-u state) x))
         (z (pdas-z state))
         (w (pdas-w state))
         (A (pdas-A state)))
    (values l u
            (matlisp:m.* w u)
            (matlisp:m.* z l)
            (matlisp:m- (sparse-m* A x)
                        (pdas-b state))
            (matlisp:m- (matlisp:m+! z
                                     (sparse-m* A (pdas-y state)
                                                :transpose t))
                        (matlisp:m+ w (pdas-c state))))))

(defun direction (state)
  (declare (type pdas-state state))
  (let ((z (pdas-z state))
        (w (pdas-w state))
        (A (pdas-A state)))
    (multiple-value-bind (l u
                          wu zl Ax-b dual)
        (violation state)
      (solve-kkt-newton l u
                        (matlisp:copy w)
                        (matlisp:copy z)
                        A
                        wu zl Ax-b dual))))

(defun box-step (l u dx)
  (let ((max double-float-positive-infinity))
    (map nil (lambda (l dx u)
               (setf dx (- dx))
               (assert (plusp l))
               (assert (plusp u))
               (cond ((zerop dx))
                     ((< dx 0)
                      (setf max (min max (/ l (- dx)))))
                     (t
                      (setf max (min max (/ u dx))))))
         (matlisp::store l)
         (matlisp::store dx)
         (matlisp::store u))
    max))

(defun pos-step (v dv)
  (let ((max double-float-positive-infinity))
    (map nil (lambda (v dv)
               (setf dv (- dv))
               (when (< dv 0)
                 (assert (> v 0))
                 (setf max (min max
                                (- (/ v dv))))))
         (matlisp::store v)
         (matlisp::store dv))
    max))

(defun pdas-step (l u state dw dx dy dz)
  (declare (ignore dy))
  (min (box-step l u dx)
       (pos-step (pdas-w state) dw)
       (pos-step (pdas-z state) dz)))

(defun apply-step (state step dw dx dy dz)
  (flet ((dec (v dv)
           (matlisp:axpy! (- step) dv v)))
    (dec (pdas-w state) dw)
    (dec (pdas-x state) dx)
    (dec (pdas-y state) dy)
    (dec (pdas-z state) dz))
  state)

(defun clear-step (state violations dw dz)
  (when (every (lambda (x)
                 (>= x *clamp*))
               (matlisp::store (pdas-u state)))
    (setf (third violations) 0d0)
    (matlisp:fill-matrix (pdas-w state) 1d-10)
    (matlisp:fill-matrix dw 0))
  (when (every (lambda (x)
                 (<= x (- *clamp*)))
               (matlisp::store (pdas-l state)))
    (setf (fourth violations) 0d0)
    (matlisp:fill-matrix (pdas-z state) 1d-10)
    (matlisp:fill-matrix dz 0)))

(defun cholesky-ls! (A s x)
  (declare (type (alien (* cholmod-sparse)) A)
           (optimize debug))
  (assert (<= (slot A 'nrow) (slot A 'ncol)))
  (with-alien ((N (* cholmod-sparse) :local
                  (cholmod-copy-sparse A *cholmod-common*)))
    (scale-sparse! N s)
    (prog1 (let ((N^-1x (solve-sparse N x)))
             (and N^-1x
                  (sparse-m* N N^-1x :transpose t)))
      (cholmod-free-sparse (addr N) *cholmod-common*))))

(defun slack (l x u max)
  (matlisp:make-real-matrix
   (map '(simple-array double-float 1)
        (lambda (l x u)
          (let ((distance (min max
                               (- x l)
                               (- u x))))
            (assert (typep distance '(double-float (0d0))))
            distance))
        (matlisp::store l)
        (matlisp::store x)
        (matlisp::store u))))

(defun residual (state)
  (matlisp:m- (pdas-b state)
              (sparse-m* (pdas-a state)
                         (pdas-x state))))

(defun max-step (l x u g)
  (let ((max-step double-float-positive-infinity))
    (dotimes (i (matlisp:nrows x) max-step)
      (let ((x (matlisp:mref x i))
            (l (matlisp:mref l i))
            (u (matlisp:mref u i))
            (g (matlisp:mref g i)))
        (let ((step (cond ((zerop g)
                           double-float-positive-infinity)
                          ((minusp g)
                           (/ (- l x) g))
                          (t
                           (/ (- u x) g)))))
          (setf max-step (min max-step step)))))))

(defun one-repair-iteration (state &optional residual)
  (let* ((x (pdas-x state))
         (l (pdas-l state))
         (u (pdas-u state))
         (slack (slack l x u 1d4))
         (residual (or residual
                       (residual state)))
         (dg (cholesky-ls! (pdas-A state)
                           slack
                           (matlisp:copy residual)))
         (g (matlisp:m.* dg slack)))
    (let* ((gamma 0.9d0)
           (step (* gamma (min (max-step l x u g) (/ gamma))))
           (norm-g (matlisp:norm g)))
      (format t "| ~5,2g ~5,2g ~%" norm-g
              (* step norm-g))
      (setf (pdas-x state)
            (matlisp:map-matrix! (lambda (x)
                                   (max x 1d-4))
                                 (matlisp:axpy step g x))))
    (values state t)))

(defun centering-direction (l x u)
  (matlisp:make-real-matrix
   (map '(simple-array double-float 1)
        (lambda (l x u)
          (cond ((and (float-infinity-p l)
                      (float-infinity-p u))
                 0d0)
                ((< (- x l) (- u x))
                 (min 1d0 (- u x)))
                (t
                 (max -1d0 (- l x)))))
        (matlisp::store l)
        (matlisp::store x)
        (matlisp::store u))))

(defun primal-project (scale c constraints)
  (with-alien ((A (* cholmod-sparse) :local
                  (cholmod-copy-sparse constraints *cholmod-common*)))
    (scale-sparse! A scale)
    (let* ((sc (matlisp:m.*! scale (matlisp:scal -1 c)))
           (Asc (sparse-m* A sc))
           (proj (solve-sparse A Asc)))
      (prog1 (and proj
                  (sparse-m* A proj :transpose t
                                    :y sc
                                    :alpha -1d0  :beta 1d0
                                    :output sc))
        (cholmod-free-sparse (addr A) *cholmod-common*)))))

(defun one-pdas-iteration (state repair)
  (multiple-value-bind (l u
                        wu zl Ax-b dual)
      (violation state)
    (assert (every #'plusp (matlisp::store l)))
    (assert (every #'plusp (matlisp::store u)))
    (let* ((pobj (matlisp:dot (pdas-c state) (pdas-x state)))
           (dobj (+ (matlisp:dot (pdas-b state) (pdas-y state))
                   (matlisp:dot (pdas-l state) (pdas-z state))
                   (- (matlisp:dot (pdas-u state) (pdas-w state)))))
           (violations (list (matlisp:norm Ax-b :inf)
                             (matlisp:norm dual :inf)
                             (matlisp:norm wu :inf)
                             (matlisp:norm zl :inf)))
           (primal-feasible (< (first violations) 1d-2))
           (feasible (and (< (first violations) 1d-2)
                          (< (second violations) 1d-2))))
      (format t "~8,3g ~8,3g | ~{~8,3g ~}"
              pobj dobj
              violations)
      (setf ; primal-feasible t
            feasible t
       )
      (unless primal-feasible
        (one-repair-iteration state)
        (return-from one-pdas-iteration
          (values (/ (abs (- pobj dobj))
                     (max (abs pobj) (abs dobj) 1d0))
                  dobj)))
      (when repair
        (matlisp:m+! 1d-4 (pdas-w state))
        (matlisp:m+! 1d-4 (pdas-z state))
        (let* ((l (pdas-l state))
               (x (pdas-x state))
               (u (pdas-u state))
               (slack (slack l x u 1d4))
               (dx (primal-project slack (centering-direction l
                                                              x
                                                              u)
                                   (pdas-A state)))
               (dx (matlisp:m.* dx slack))
               (step (* .5d0 (max-step l x u dx))))
          (setf (pdas-x state)
                (matlisp:axpy! step dx (pdas-x state))))
        (return-from one-pdas-iteration
          (values (/ (abs (- pobj dobj))
                     (max (abs pobj) (abs dobj) 1d0))
                  dobj)))
      (multiple-value-bind (dw dx dy dz)
          (solve-kkt-newton
           (matlisp:m+ 0 l) (matlisp:m+ 0 u)
           (matlisp:m+ 0 (pdas-w state))
           (matlisp:m+ 0 (pdas-z state))
           (pdas-A state)
           wu zl Ax-b dual)
        #+nil(clear-step state violations dw dz)
        (let ((step (pdas-step l u state dw dx dy dz)))
          (format t "| ~8,3g ~%" step)
          (apply-step state (the (double-float (0d0) 1d0)
                                 (min 1d0 (* .9d0 step)))
                      dw dx dy dz)
          (values (/ (abs (- pobj dobj))
                     (max (abs pobj) (abs dobj) 1d0))
                  dobj
                  step))))))

(defun pdas (state &optional (max-iter most-positive-fixnum))
  (with-cholmod ()
    (unwind-protect (loop with repair = nil
                          for i upfrom 1 upto max-iter do
                          (format t "~&~3d:~10t" i)
                          (multiple-value-bind (viol obj step)
                              (one-pdas-iteration state
                                                  repair)
                            (setf repair (and step (< step 1d-6)))
                            (when (< viol 1d-4)
                              (return (values obj viol)))))
      (free-pdas-a state))))
