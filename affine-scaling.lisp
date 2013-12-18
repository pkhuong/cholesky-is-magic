(defstruct (affine-scaling-state
            (:conc-name #:affine-))
  nvars ncons
  x
  c
  triplets
  %A %A-copy
  b
  l u
  (sparse-state (make-solve-sparse-state) :read-only t))

(defun affine-A (affine)
  (or (affine-%a affine)
      (let ((A (make-sparse-from-triplet-vector
                                (affine-ncons affine)
                                (affine-nvars affine)
                                (affine-triplets affine))))
        (assert (not (affine-%a-copy affine)))
        (setf (affine-%a-copy affine) (cholmod-copy-sparse
                                       A
                                       *cholmod-common*)
              (affine-%a affine) A))))

(define-alien-routine memcpy (* t)
  (dst (* t))
  (src (* t))
  (len size-t))

(defun affine-A-copy (affine)
  (let ((A (affine-A affine))
        (copy (affine-%a-copy affine)))
    (declare (type (alien (* cholmod-sparse)) A copy))
    (memcpy (slot copy 'x) (slot a 'x)
            (* 8 (slot a 'nzmax)))
    copy))

(defun free-sparse (sparse)
  (with-alien ((sparse (* cholmod-sparse)
                       :local sparse))
    (assert (/= (cholmod-free-sparse (addr sparse)
                                     *cholmod-common*)
                0))))

(defun free-affine-A (affine)
  (let ((sparse (shiftf (affine-%a affine) nil))
        (copy   (shiftf (affine-%a-copy affine) nil)))
    (when sparse
      (free-sparse sparse))
    (when copy
      (free-sparse copy))))

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

(defun solve-symmetric! (A x &optional sparse-state)
  ;; A\x, destructive on both
  ;;
  ;; Cholesky *should* be Good Enough.
  (solve-sparse A x sparse-state))

(defun project (scale c constraints sparse-state)
  ;; min || x + [scale]c ||_2
  ;;  s.t.  [constraints][scale]x = 0
  ;;
  ;; Let B = [constraints][scale], N = BB', c = -[scale]c
  ;;
  ;; x = c - B'N^-1 B c
  (let* ((sc (matlisp:m.*! scale
                           (matlisp:scal -1 c)))
         (AD (scale-sparse! constraints scale))
         (AD2c (sparse-m* AD sc))
         (N^-1AD2c (solve-symmetric! AD AD2c
                                     sparse-state)))
    ;; Dc - (AD)t N^-1 AD Dc
    (and N^-1AD2c
         (values (sparse-m* AD N^-1AD2c :transpose t :y sc
                                        :alpha -1d0 :beta 1d0
                                        :output sc)
                 N^-1AD2c))))

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

(defun one-affine-scaling-iteration (state &key centering)
  (declare (type affine-scaling-state state))
  (let* ((x (affine-x state))
         (l (affine-l state))
         (u (affine-u state))
         (slack (slack l x u *max-slack*)))
    (multiple-value-bind (dg ye)
        (project slack (if centering
                           (centering-direction l x u)
                           (affine-c state))
                 (affine-A-copy state)
                 (affine-sparse-state state))
      (declare (ignore ye))
      (unless dg
        (format t " singular ")
        (return-from one-affine-scaling-iteration
          (values state nil)))
      (let* ((g (matlisp:m.* dg slack))
             (step (* *gamma* (max-step l x u g)))
             (norm-g (matlisp:norm g))
             (norm-dg (matlisp:norm dg))
             (descent (matlisp:dot g (affine-c state))))
        (when (> step 1d10)
          (error "Unbounded problem"))
        (format t "~12,5g ~12,5g"
                (* step norm-g)
                norm-dg)
        (unless centering
          (when (or (< norm-dg (min 1d-6
                                    (* 1d-8 (matlisp:nrows x))))
                    (plusp descent))
            (when (plusp descent)
              (format t " Not a descent direction "))
            (return-from one-affine-scaling-iteration
              (values state nil)))
          (when (or (< (* step norm-g) 1d-6)
                    (plusp descent))
            (format t " ... ")
            (return-from one-affine-scaling-iteration
              (one-affine-scaling-iteration state :centering t))))
        (setf (affine-x state)
              (matlisp:axpy! step g x))))
    (values state t)))

(defun residual (state)
  (declare (type affine-scaling-state state))
  (matlisp:m- (affine-b state)
              (sparse-m* (affine-a state)
                         (affine-x state))))

(defun cholesky-ls! (A x sparse-state)
  (declare (type (alien (* cholmod-sparse)) A)
           (optimize debug))
  (assert (<= (slot A 'nrow) (slot A 'ncol)))
  (let ((N^-1x (solve-symmetric! A x sparse-state)))
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
         (constraints (scale-sparse! (affine-A-copy state) slack))
         (dg (cholesky-ls! constraints residual
                           (affine-sparse-state state)))
         (g (matlisp:m.* dg slack)))
    (let ((step (* *gamma* (min (max-step l x u g) (/ *gamma*))))
          (norm-g (matlisp:norm g)))
      (format t "~12,5g " norm-g)
      (setf (affine-x state)
            (matlisp:axpy step g x)))
    (values state t)))

(defun one-iteration (state &optional centering)
  (declare (type affine-scaling-state state))
  (let* ((residual (residual state))
         (norm (matlisp:norm residual)))
    (cond ((> norm (* 1d-6 (matlisp:nrows residual)))
           (prog2 (format t "Repair:   ~12,5g ... " norm)
               (multiple-value-list (one-repair-iteration state residual))
             (format t "~12,5g~%" (matlisp:norm (residual state)))))
          (t
           (prog2 (format t "~A ~12,5g ... "
                          (if centering
                              "Recenter:"
                              "Optimize:")
                          (matlisp:dot (affine-x state) (affine-c state)))
               (multiple-value-list (one-affine-scaling-iteration
                                     state
                                     :centering centering))
             (format t "~12,5g~%" (matlisp:dot (affine-x state)
                                             (affine-c state))))))))

(defun affine-scaling (state)
  (declare (type affine-scaling-state state))
  (with-cholmod ()
    (unwind-protect
         (let ((common *cholmod-common*))
           (setf (solve-sparse-state-factor (affine-sparse-state state))
                 (cholmod-analyze (affine-A state) common))
           (assert (zerop (cholmod-get-status common)))
           (format t "~
AA':    nnz: ~12,5g flops: ~12,5g
Factor: nnz: ~12,5g flops: ~12,5g~%"
                   (cholmod-get-anz common)
                   (cholmod-get-aatfl common)
                   (cholmod-get-lnz common)
                   (cholmod-get-fl common))
           (loop for i upfrom 0 do
             (format t "~4d: " i)
             (destructuring-bind (state continue)
                 (one-iteration state (zerop (mod (1+ i) 16)))
               (let* ((residual (residual state))
                      (norm (matlisp:norm residual)))
                 (unless (or continue
                             (> norm (* 1d-6 (matlisp:nrows residual))))
                   (return (values (matlisp:dot (affine-x state)
                                                (affine-c state))
                                   (affine-x state)
                                   residual)))))))
      (free-affine-a state)
      (free-sparse-state (affine-sparse-state state))
      (cholmod-free-work *cholmod-common*)
      (format t "Remaining space: ~A ~A ~%"
              (cholmod-get-malloc-count *cholmod-common*)
              (cholmod-get-memory-inuse *cholmod-common*)))))
