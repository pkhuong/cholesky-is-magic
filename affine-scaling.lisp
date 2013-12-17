(defstruct (affine-scaling-state
            (:conc-name #:affine-))
  x
  c
  A b
  l u)

(defun make-affine-state (sf)
  (declare (type standard-form sf))
  (let ((A (matlisp:make-real-matrix-dim (sf-ncons sf)
                                         (sf-nvars sf)))
        (x (matlisp:make-real-matrix-dim (sf-nvars sf) 1))
        (l (matlisp:make-real-matrix (sf-l sf)))
        (u (matlisp:make-real-matrix (sf-u sf))))
    (map nil (lambda (triplet)
               (declare (type triplet triplet))
               (setf (matlisp:matrix-ref A
                                         (triplet-row triplet)
                                         (triplet-col triplet))
                     (triplet-value triplet)))
         (sf-A sf))
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
     :x x
     :c (let ((v (matlisp:make-real-matrix-dim (sf-nvars sf) 1)))
          (loop for (i . x) across (sf-c sf)
                do (setf (matlisp:matrix-ref v i) x))
          v)
     :A A
     :b (matlisp:make-real-matrix (sf-b sf))
     :l l
     :u u)))

(defstruct (project-work-space
            (:conc-name #:proj-))
  A
  scale
  B
  BBt
  c
  d
  lwork
  work)

(defvar *work-space*)

(defun allocate-work-space (n p)
  (or *work-space*
      (setf *work-space*
            (let ((lwork (+ n p (* (max p n) 10))))
              (make-project-work-space
               :A (matlisp:make-real-matrix-dim n n)
               :scale (matlisp:make-real-matrix-dim n n)
               :B (matlisp:make-real-matrix-dim p n)
               :BBt (matlisp:make-real-matrix-dim p p)
               :c (matlisp:make-real-matrix-dim n 1)
               :d (matlisp:make-real-matrix-dim p 1)
               :lwork lwork
               :work (matlisp::allocate-real-store lwork))))))

(defun scale-constraints (constraints scale)
  (let* ((n (matlisp:ncols constraints))
         (p (matlisp:nrows constraints))
         (proj (allocate-work-space n p))
         (scale (setf (matlisp:diag (proj-scale proj)) scale)))
    (matlisp:gemm! 1 constraints scale
                   0 (proj-B proj))))

(defun solve-symmetric! (A x &aux (n (matlisp:nrows A)))
  ;; A\x, destructive on both
  ;;
  ;; Cholesky *should* be Good Enough.
  (multiple-value-bind (factor info)
      (lapack:dpotrf "L" n
                     (matlisp::store A) (max 1 n)
                     0)
    (unless (zerop info)
      (format t " (DPOTRF: ~A) " info)
      (return-from solve-symmetric!))
    (multiple-value-bind (solution info)
        (lapack:dpotrs "L" n 1
                       factor (max 1 n)
                       (matlisp::store x) (max 1 n)
                       0)
      (assert (zerop info))
      (setf (matlisp::store x) solution)
      x)))

(defun project (scale c constraints)
  ;; min || x + [scale]c ||_2
  ;;  s.t.  [constraints][scale]x = 0
  ;;
  ;; Let B = [constraints][scale], N = BB', c = -[scale]c
  ;;
  ;; x = c - B'N^-1 B c
  (let* ((m (matlisp:nrows c))
         (n m)
         (p (matlisp:nrows constraints))
         (proj (allocate-work-space n p))
         (c (matlisp:m.*! scale
                          (matlisp:scal -1 c)))
         (B (scale-constraints constraints scale))
         (N (matlisp:gemm! 1 B B 0 (proj-BBt proj) :nt))
         (Bc (matlisp:m* B c))
         (N^-1Bc (or (solve-symmetric! N bc)
                     (matlisp:gelsy! (matlisp:gemm! 1 B B 0 (proj-BBt proj)
                                                    :nt)
                                     (matlisp:m* B c)
                                     1d-10))))
    ;; c - Bt N^-1 B c
    (matlisp:gemm! -1 B N^-1Bc 1 c :tn)))

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
              (matlisp:m* (affine-a state)
                          (affine-x state))))

(defun cholesky-ls! (A x)
  (declare (optimize debug))
  (assert (<= (matlisp:nrows A) (matlisp:ncols A)))
  (let* ((n (matlisp:ncols A))
         (out (matlisp:make-real-matrix-dim n 1))
         (proj (allocate-work-space n (matlisp:nrows A)))
         (N (matlisp:gemm! 1 A A 0 (proj-BBt proj) :nt))
         (N^-1x (solve-symmetric! N x)))
    (and N^-1x
         (matlisp:gemm! 1 A N^-1x
                        0 out
                        :tn))))

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
         (constraints (affine-A state))
         (dg (or (cholesky-ls! (scale-constraints constraints slack)
                               (matlisp:copy residual))
                 (matlisp:gelsy! (scale-constraints constraints slack)
                                 residual 1d-12)))
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
                                residual)))))))
