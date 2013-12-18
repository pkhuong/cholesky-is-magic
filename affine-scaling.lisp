(defstruct (affine-scaling-state
            (:conc-name #:affine-))
  nvars ncons
  x
  c
  triplets
  %A b
  l u)

(defun affine-A (affine)
  (or (affine-%a affine)
      (setf (affine-%a affine)
            (make-sparse-from-triplet-vector
             (affine-ncons affine)
             (affine-nvars affine)
             (affine-triplets affine)))))

(defun free-affine-A (affine)
  (let ((sparse (shiftf (affine-%a affine) nil)))
    (when sparse
      (with-alien ((sparse (* cholmod-sparse)
                           :local sparse))
        (assert (/= (cholmod-free-sparse (addr sparse)
                                         *cholmod-common*)
                    0))))))

(defun make-affine-state (sf)
  (declare (type standard-form sf))
  (let ((x (matlisp:make-real-matrix-dim (sf-nvars sf) 1))
        (l (matlisp:make-real-matrix (sf-l sf)))
        (u (matlisp:make-real-matrix (sf-u sf))))
    (dotimes (i (sf-nvars sf))
      (when (< (- (matlisp:matrix-ref u i)
                  (matlisp:matrix-ref l i))
               1d-6)
        (decf (matlisp:matrix-ref l i) 5d-7)
        (incf (matlisp:matrix-ref u i) 5d7))
      (let* ((l (matlisp:matrix-ref l i))
             (u (matlisp:matrix-ref u i))
             (delta (- u l)))
        (setf (matlisp:matrix-ref x i)
              (cond ((and (< l -1d10) (> u 1d10))
                       
                     0d0)
                    ((< l -1d10)
                     (- u (min (/ delta 2)
                               (1+ (* (abs u) .1d0)))))
                    ((> u 1d10)
                     (+ l (min (/ delta 2)
                               (1+ (* (abs l) 1d0)))))
                    (t
                     (/ (+ l u) 2))))))
    (make-affine-scaling-state
     :nvars (sf-nvars sf)
     :ncons (sf-ncons sf)
     :x x
     :c (let ((v (matlisp:make-real-matrix-dim (sf-nvars sf) 1)))
          (loop for (i . x) across (sf-c sf)
                do (setf (matlisp:matrix-ref v i) x))
          v)
     :triplets (sf-A sf)
     :%A nil
     :b (matlisp:make-real-matrix (sf-b sf))
     :l l
     :u u)))

(defstruct (project-work-space
            (:conc-name #:proj-))
  scale
  B
  BBt
  c
  d)

(defvar *work-space*)

(defun allocate-work-space (n p)
  (or *work-space*
      (setf *work-space*
            (make-project-work-space
             :scale (matlisp:make-real-matrix-dim n n)
             :B (matlisp:make-real-matrix-dim p n)
             :BBt (matlisp:make-real-matrix-dim p p)
             :c (matlisp:make-real-matrix-dim n 1)
             :d (matlisp:make-real-matrix-dim p 1)))))

(defun scale-constraints (constraints scale)
  (let* ((n (matlisp:ncols constraints))
         (p (matlisp:nrows constraints))
         (proj (allocate-work-space n p))
         (scale (setf (matlisp:diag (proj-scale proj)) scale)))
    (matlisp:gemm! 1 constraints scale
                   0 (proj-B proj))))

(defun solve-symmetric! (A x)
  ;; A\x, destructive on both
  ;;
  ;; Cholesky *should* be Good Enough.
  (solve-sparse A x))

(defun project (scale c constraints)
  ;; min || x + [scale]c ||_2
  ;;  s.t.  [constraints][scale]x = 0
  ;;
  ;; Let B = [constraints][scale], N = BB', c = -[scale]c
  ;;
  ;; x = c - B'N^-1 B c
  (let* ((c (matlisp:m.*! scale
                          (matlisp:scal -1 c)))
         (B (scale-sparse constraints scale))
         (Bc (sparse-m* B c))
         (N^-1Bc (solve-symmetric! B bc)))
    ;; c - Bt N^-1 B c
    (prog1 (sparse-m* B N^-1Bc :transpose t :y c
                               :alpha -1d0 :beta 1d0)
      (with-alien ((B (* cholmod-sparse) :local B))
        (cholmod-free-sparse (addr B) *cholmod-common*)))))

(defvar *max-slack* 1d8)

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

(defvar *gamma* .9d0)

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

(defun one-affine-scaling-iteration (state)
  (declare (type affine-scaling-state state))
  (let* ((x (affine-x state))
         (l (affine-l state))
         (u (affine-u state))
         (slack (slack l x u *max-slack*))
         (dg (project slack
                      (affine-c state)
                      (affine-A state)))
         (g (matlisp:m.* dg slack)))
    (let ((step (* *gamma* (max-step l x u g)))
          (norm-g (matlisp:norm g)))
      (when (> step 1d10)
        (error "Unbounded problem"))
      (format t "~12,5g " norm-g)
      (when (or (< (min (* step norm-g) norm-g)
                   (min (* 1d-10 (matlisp:nrows g))
                        1d-6))
                (plusp (matlisp:dot g (affine-c state))))
        (return-from one-affine-scaling-iteration
          (values state nil)))
      (setf (affine-x state)
            (matlisp:axpy! step g x)))
    (values state t)))

(defun residual (state)
  (declare (type affine-scaling-state state))
  (matlisp:m- (affine-b state)
              (sparse-m* (affine-a state)
                         (affine-x state))))

(defun cholesky-ls! (A x)
  (declare (type (alien (* cholmod-sparse)) A)
           (optimize debug))
  (assert (<= (slot A 'nrow) (slot A 'ncol)))
  (let ((N^-1x (solve-symmetric! A x)))
    (and N^-1x
         (sparse-m* A N^-1x :transpose t))))

;; Least-square solve for
;; Asx - r = 0
;;  min ||x||
(defun one-repair-iteration (state &optional residual)
  (declare (type affine-scaling-state state))
  (let* ((x (affine-x state))
         (l (affine-l state))
         (u (affine-u state))
         (slack (slack l x u (sqrt *max-slack*)))
         (residual (or residual
                        (residual state)))
         (constraints (scale-sparse (affine-A state) slack))
         (dg (cholesky-ls! constraints
                           (matlisp:copy residual)))
         (g (matlisp:m.* dg slack)))
    (with-alien ((constraints (* cholmod-sparse)
                              :local constraints))
      (cholmod-free-sparse (addr constraints) *cholmod-common*))
    (let ((step (* *gamma* (min (max-step l x u g) (/ *gamma*))))
          (norm-g (matlisp:norm g)))
      (format t "~12,5g " norm-g)
      (setf (affine-x state)
            (matlisp:axpy step g x)))
    (values state t)))

(defun one-iteration (state)
  (declare (type affine-scaling-state state))
  (let* ((residual (residual state))
         (norm (matlisp:norm residual)))
    (cond ((> norm (* 1d-6 (matlisp:nrows residual)))
           (prog2 (format t "Repair:   ~12,5g ... " norm)
               (multiple-value-list (one-repair-iteration state residual))
             (format t "~12,5g~%" (matlisp:norm (residual state)))))
          (t
           (prog2 (format t "Optimise: ~12,5g ... "
                          (matlisp:dot (affine-x state) (affine-c state)))
               (multiple-value-list (one-affine-scaling-iteration state))
             (format t "~12,5g~%" (matlisp:dot (affine-x state)
                                             (affine-c state))))))))

(defun affine-scaling (state)
  (declare (type affine-scaling-state state))
  (with-cholmod ()
    (unwind-protect
         (loop with *work-space* = nil
               for i upfrom 0 do
                 (format t "~4d: " i)
                 (destructuring-bind (state continue)
                     (one-iteration state)
                   (let* ((residual (residual state))
                          (norm (matlisp:norm residual)))
                     (unless (or continue
                                 (> norm (* 1d-6 (matlisp:nrows residual))))
                       (return (values (matlisp:dot (affine-x state)
                                                    (affine-c state))
                                       (affine-x state)
                                       residual))))))
      (free-affine-a state)
      (cholmod-free-work *cholmod-common*)
      (format t "Remaining space: ~A ~A ~%"
              (cholmod-get-malloc-count *cholmod-common*)
              (cholmod-get-memory-inuse *cholmod-common*)))))
