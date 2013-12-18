(defstruct (affine-scaling-state
            (:conc-name #:affine-))
  nvars ncons
  x
  c
  triplets
  %A %A-copy
  b
  l u
  (sparse-state (make-solve-sparse-state)))

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

(defun scale-constraints (constraints scale)
  (let* ((n (matlisp:ncols constraints))
         (p (matlisp:nrows constraints))
         (proj (allocate-work-space n p))
         (scale (setf (matlisp:diag (proj-scale proj)) scale)))
    (matlisp:gemm! 1 constraints scale
                   0 (proj-B proj))))

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
  (let* ((c (matlisp:m.*! scale
                          (matlisp:scal -1 c)))
         (B (scale-sparse! constraints scale))
         (Bc (sparse-m* B c))
         (N^-1Bc (solve-symmetric! B bc sparse-state)))
    ;; c - Bt N^-1 B c
    (sparse-m* B N^-1Bc :transpose t :y c
                        :alpha -1d0 :beta 1d0)))

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
                      (affine-A-copy state)
                      (affine-sparse-state state)))
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
      (free-sparse-state (affine-sparse-state state))
      (cholmod-free-work *cholmod-common*)
      (format t "Remaining space: ~A ~A ~%"
              (cholmod-get-malloc-count *cholmod-common*)
              (cholmod-get-memory-inuse *cholmod-common*)))))
