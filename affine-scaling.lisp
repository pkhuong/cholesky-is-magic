(defstruct (affine-scaling-state
            (:conc-name #:affine-))
  x
  c
  A b
  l u)

(defun make-affine-state (sf)
  (declare (type standard-form sf))
  (let ((A (matlisp:make-real-matrix-dim (sf-ncons sf)
                                         (sf-nvars sf))))
    (map nil (lambda (triplet)
               (declare (type triplet triplet))
               (setf (matlisp:matrix-ref A
                                         (triplet-row triplet)
                                         (triplet-col triplet))
                     (triplet-value triplet)))
         (sf-A sf))
    (make-affine-scaling-state
     :x (matlisp:make-real-matrix
         (map '(simple-array double-float 1)
              (lambda (l u)
                (cond ((< l -1d10)
                       (1- u))
                      ((> u 1d10)
                       (1+ l))
                      (t
                       (/ (+ l u) 2))))
              (sf-l sf)
              (sf-u sf)))
     :c (let ((v (matlisp:make-real-matrix-dim (sf-nvars sf) 1)))
          (loop for (i . x) across (sf-c sf)
                do (setf (matlisp:matrix-ref v i) x))
          v)
     :A A
     :b (matlisp:make-real-matrix (sf-b sf))
     :l (matlisp:make-real-matrix (sf-l sf))
     :u (matlisp:make-real-matrix (sf-u sf)))))

(defun project (scale c constraints &optional rhs)
  ;; min || Ix - [scale]c ||_2
  ;;  s.t.  [constraints][scale]x = 0
  ;;
  ;; In DGGLSE:
  ;;   A = I
  ;;   B = A[scale]
  (let* ((m (matlisp:nrows c))
         (n m)
         (p (matlisp:nrows constraints))
         (A (matlisp:scal! -1 (matlisp:eye n)))
         (scale (let ((m (matlisp:make-real-matrix-dim n m)))
                  (setf (matlisp:diag m) scale)
                  m))
         (B (matlisp:m* constraints scale))
         (c (matlisp:m* scale c))
         (d (if rhs
                (matlisp:copy rhs)
                (matlisp:make-real-matrix-dim p 1)))
         (lwork (+ n p (* (max m p n) 10))))
    (multiple-value-bind (A B c d x work info)
        (lapack:dgglse m n p
                       (matlisp::store A) (max 1 m)
                       (matlisp::store B) (max 1 p)
                       (matlisp::store c) (matlisp::store d)
                       (matlisp::allocate-real-store m)
                       (matlisp::allocate-real-store lwork) lwork
                       0)
      (declare (ignore A B c d))
      (assert (= info 0))
      (let ((x (make-instance 'matlisp:real-matrix
                              :nrows m
                              :ncols 1
                              :store x)))
        (values x (ceiling (aref work 0)))))))

(defvar *max-slack* 1d10)

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

(defvar *gamma* 98d-2)

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
    ;; better be a descent direction
    (assert (minusp (matlisp:dot g (affine-c state))))
    ;; and better satisfy the equality constraint
    (assert (< (matlisp:norm (matlisp:m* (affine-A state) g))
               1d-4))
    (let ((step (* *gamma* (max-step l x u g)))
          (norm-g (matlisp:norm g)))
      (when (> step 1d10)
        (error "Unbounded problem"))
      (format t "~12,5g " norm-g)
      (when (< (min (* step norm-g) norm-g) 1d-6)
        (return-from one-affine-scaling-iteration
          (values state nil)))
      (setf (affine-x state)
            (matlisp:axpy! step g x)))
    (values state t)))

(defun residual (state)
  (declare (type affine-scaling-state state))
  (matlisp:m- (affine-b state)
              (matlisp:m* (affine-a state)
                          (affine-x state))))

;; A(w [dg]) = r

(defun one-repair-iteration (state &optional residual)
  (declare (type affine-scaling-state state))
  (let* ((x (affine-x state))
         (l (affine-l state))
         (u (affine-u state))
         (slack (slack l x u (sqrt *max-slack*)))
         (residual (or residual
                        (residual state)))
         (dg (project slack
                      (matlisp:make-real-matrix-dim (matlisp:nrows x) 1)
                      (affine-A state)
                      residual))
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
    (cond ((> norm 5d-7)
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
  (loop for i upfrom 0 do
    (format t "~4d: " i)
    (destructuring-bind (state continue)
        (one-iteration state)
      (unless continue
        (return (values (matlisp:dot (affine-x state)
                                     (affine-c state))
                        (affine-x state)
                        (matlisp:norm (residual state))))))))
