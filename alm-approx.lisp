;;; Accelerated parallel proximal coordinate descent
;;;
;;; For now, assume that we perform full-vector descent, without
;;; any random sampling.

#||
Minimise f(x) = sum_j f_j(x)

subject to l <= x <= u

||#
(deftype index ()
  `(mod ,most-positive-fixnum))

(deftype ivec (&optional length)
  `(simple-array index (,length)))

(deftype dvec (&optional length)
  `(simple-array double-float (,length)))

;; linear obj -> Lipschitz constant = 0
(defstruct (linear
            (:constructor %make-linear))
  (indices (error "Missing arg") :type ivec)
  (coefs (error "Missing arg") :type dvec))

(defun make-linear (pairs)
  (let ((pairs (remove 0 pairs :key #'cdr :test #'=)))
    (%make-linear
     :indices (map 'ivec #'car pairs)
     :coefs (map 'dvec (lambda (x)
                         (float (cdr x) 1d0))
                 pairs))))

;; 1/2 |ax-b|^2
(defstruct (quadratic
            (:constructor %make-quadratic)
            (:conc-name #:quad-))
  (scale 1d0 :type double-float)
  (indices (error "Missing arg") :type ivec)
  (coefs (error "Missing arg") :type dvec)
  (rhs (error "Missing arg") :type double-float)
  (beta (error "Missing arg") :type double-float))

(defun make-quadratic (pairs &key rhs scale)
  (let ((pairs (remove 0 pairs :key #'cdr :test #'=)))
    (%make-quadratic
     :scale (or scale 1d0)
     :indices (map 'ivec #'car pairs)
     :coefs (map 'dvec (lambda (x)
                         (float (cdr x) 1d0))
                 pairs)
     :rhs (float (or rhs 0d0) 1d0)
     ;; specialised version for tau = n
     :beta (float (length pairs) 1d0))))

(defun 2norm (values &key initial-value)
  (sqrt (reduce #'+ values :key (lambda (x)
                                  (expt x 2))
                           :initial-value (if initial-value
                                              (expt initial-value 2)
                                              0d0))))

(defun inf-norm (values)
  (reduce #'max values :key #'abs))

(defun scale-quadratic (q)
  (let ((norm (2norm (quad-coefs q)
                     :initial-value (quad-rhs q))))
    (when (> norm 1d-6)
      (let ((scale (/ norm)))
        (setf (quad-rhs q) (* (quad-rhs q) scale))
        (map-into (quad-coefs q)
                  (lambda (c)
                    (* scale c))
                  (quad-coefs q)))))
  q)

(declaim (maybe-inline violation))
(defun violation (quadratic x &optional (scale t))
  "Compute ax-b"
  (let ((acc (- (quad-rhs quadratic))))
    (declare (type double-float acc))
    (map nil (lambda (i c)
               (incf acc (* c (aref x i))))
         (quad-indices quadratic)
         (quad-coefs quadratic))
    (* acc (if scale (quad-scale quadratic) 1d0))))


;; Global Lipschitz-esque value
(defgeneric accumulate-nu (constraint nu)
  (:method ((c linear) nu)
    ;; Lipschitz constant = 0: we're elementwise
    nu)
  (:method ((q quadratic) nu)
    (let ((beta (quad-beta q))
          (scale^2 (expt (quad-scale q) 2)))
      (map nil (lambda (i c)
                 (incf (aref nu i) (* beta scale^2 c c)))
           (quad-indices q)
           (quad-coefs q)))
    nu))

(defstruct (approx-state
            (:conc-name #:approx-)
            (:constructor %make-approx-state))
  (orig-vars (error "Missing arg") :type index)
  (orig-cons (error "Missing arg") :type index)
  (nvars (error "Missing arg") :type index)
  (constraints (error "Missing arg") :type simple-vector)
  (nu (error "Missing arg") :type dvec)
  (l (error "Missing arg") :type dvec)
  (u (error "Missing arg") :type dvec)
  (z0 0d0 :type double-float)
  (weight (error "Missing arg") :type double-float)
  (lambda (error "Missing arg") :type dvec)
  )

(defun make-approx-state (orig-vars orig-cons
                          constraints nvars l u
                          weight lambda z0)
  (let* ((constraints (coerce constraints 'simple-vector))
         (nu (make-array nvars :element-type 'double-float
                               :initial-element 0d0)))
    (map nil (lambda (c)
               (accumulate-nu c nu))
         constraints)
    (%make-approx-state :orig-vars orig-vars :orig-cons orig-cons
                        :nvars nvars
                        :constraints constraints
                        :nu nu
                        :l l
                        :u u
                        :weight weight
                        :lambda lambda
                        :z0 z0)))

(defun dual-value (approx x)
  (+ (approx-z0 approx)
     (%value-&-gradient (aref (approx-constraints approx)
                              (approx-orig-cons approx))
                        x)))

(defun primal-value (approx x)
  (+ (approx-z0 approx)
     (value-&-gradient approx x)))

(defgeneric %value-&-gradient (constraint x &optional g)
  (:method ((c linear) x &optional g)
    (declare (type dvec x)
             (type (or null dvec) g)
             (optimize speed))
    (let ((acc 0d0))
      (declare (type double-float acc))
      (map nil (lambda (i v)
                 (incf acc (* v (aref x i)))
                 (when g
                   (incf (aref g i) v)))
           (linear-indices c)
           (linear-coefs c))
      acc))
  (:method ((q quadratic) x &optional g)
    (declare (type dvec x)
             (type (or null dvec) g)
             (optimize speed)
             (inline violation))
    (let ((violation (violation q x)) ;; ax-b
          (scale (quad-scale q)))
      (when g
        (map nil (lambda (i c)
                   (incf (aref g i) (* c scale violation)))
             (quad-indices q)
             (quad-coefs q)))
      (* .5d0 (expt violation 2)))))

(defun value-&-gradient (state x)
  (declare (optimize speed))
  (let ((g (make-array (approx-nvars state)
                       :element-type 'double-float
                       :initial-element 0d0))
        (z 0d0)
        (max 0d0))
    (declare (type double-float z max))
    (map nil (lambda (c)
               (let ((v (%value-&-gradient c x g)))
                 (declare (type double-float v))
                 (incf z v)
                 (unless (linear-p c)
                   (setf max (max max (abs v))))))
         (approx-constraints state))
    (values z
            g
            max)))

(declaim (maybe-inline solve-coordinate))
(defun solve-coordinate (z nu theta g l u)
  "Min_x gx + (n theta nu)/2tau (x-z)^2, l <= x <= u,
   Tau = n"
  (declare (type double-float z nu theta g l u)
           (optimize speed))
  (let ((step (* theta nu)))
    (cond ((zerop step)
           (cond ((minusp g) u)
                 ((zerop g) z)
                 (t l)))
          (t
           (let ((best (- z (* 0.95d0 (/ g step)))))
             (cond ((< best l)
                    l)
                   ((> best u)
                    u)
                   (t best)))))))

(defun approx-descent (state theta y z)
  (declare (type double-float theta)
           (type dvec y z)
           (optimize speed)
           (inline solve-coordinate))
  (multiple-value-bind (value g max)
      (value-&-gradient state y)
    (declare (type dvec g))
    (values
     (let ((x (make-array (length z) :element-type 'double-float))
           (nu (approx-nu state))
           (l (approx-l state))
           (u (approx-u state)))
       (dotimes (i (length x) x)
         (let ((z (aref z i))
               (nu (aref nu i))
               (g (aref g i))
               (l (aref l i))
               (u (aref u i)))
           (setf (aref x i)
                 (solve-coordinate z nu theta g l u)))))
     g
     value
     max)))

(defun approx-iteration (state theta x z)
  (declare (type double-float theta)
           (type dvec x z)
           (optimize speed))
  (let* ((y (map-into (make-array (length x) :element-type 'double-float)
                      (lambda (x z)
                        (+ (* (- 1d0 theta) x)
                           (* theta z)))
                      x z)))
    (multiple-value-bind (zp g)
        (approx-descent state theta y z)
      (declare (type dvec zp g))
      (values (map-into (make-array (length y) :element-type 'double-float)
                        (lambda (y z zp)
                          (+ y (* theta (- zp z))))
                        y z zp)
              zp
              (* .5d0 (- (sqrt (let ((square (* theta theta)))
                                 (* (+ 4 square) square))
                               #+nil (+ (* theta 4)
                                        (* 4 (expt theta 2))))
                         (expt theta 2)))
              g))))

(defun project-gradient (state x g)
  (declare (type dvec x g)
           (optimize speed))
  (let ((p (make-array (length x) :element-type 'double-float)))
    (loop for i upfrom 0
          for l of-type double-float across (approx-l state)
          for x of-type double-float across x
          for g of-type double-float across g
          for u of-type double-float across (approx-u state)
          do (setf (aref p i)
                   (let ((xp (- x g)))
                     (when (< xp l)
                       (setf xp l))
                     (when (> xp u)
                       (setf xp u))
                     (- x xp)))
          finally (return p))))

(defun dot-diff (g x xp)
  (declare (type dvec g x xp)
           (optimize speed))
  (let ((acc 0d0))
    (declare (type double-float acc))
    (map nil (lambda (g x xp)
               (incf acc (* g (- xp x))))
         g x xp)
    acc))

(defun project (state x)
  (declare (type dvec x)
           (optimize speed))
  (let ((p (make-array (length x) :element-type 'double-float))
        (l (approx-l state))
        (u (approx-u state)))
    (dotimes (i (length x) p)
      (setf (aref p i)
            (let ((x (aref x i))
                  (l (aref l i))
                  (u (aref u i)))
              (min u (max l x)))))))

(defvar *approx-iterations* nil)

(defun approx (state n &optional x (accuracy 1d-5))
  (let* ((x (project state
                     (or x
                         (make-array (approx-nvars state)
                                     :element-type 'double-float
                                     :initial-element 0d0))))
         (z (copy-seq x))
         zp gy
         (theta 1d0))
    (dotimes (i n (assert nil))
      (setf (values x zp theta gy)
            (approx-iteration state theta x z))
      (multiple-value-bind (value g max)
          (value-&-gradient state zp)
        (cond ((plusp (dot-diff g z zp))
               #+nil(format t "~&~10tRestart ~d~%" (1+ i))
               (setf x z
                  theta 1d0))
              (t
               (setf z zp)))
        (let* ((pg (2norm (project-gradient state z g)))
               (done (or (and (> i 10)
                              (or (< pg accuracy)
                                  #+nil(< pg (* 1d-6 (1+ (abs value))))))
                         (= i (1- n)))))
          (when (or (zerop i)
                    (= i (1- n))
                    (zerop (mod (1+ i) 10000))
                    done)
            (format t "~&~10t~6d: ~12,5g ~12,5g ~12,5g ~12,5g ~12,5g~%"
                    (1+ i)
                    (2norm g)
                    pg
                    max
                    (+ value (approx-z0 state))
                    (dual-value state zp)))
          (when done
            (when *approx-iterations*
              (incf *approx-iterations* (1+ i)))
            (return (values z pg))))))))

(defun dot (x y)
  (let ((acc 0d0))
    (map nil (lambda (x y)
               (incf acc (* x y)))
         x y)
    acc))

(defun make-alm-subproblem (sf lambda weight)
  (let* ((nvars (sf-nvars sf))
         (ncons (sf-ncons sf))
         (constraints (make-array (1+ ncons) :initial-element nil))
         (l (make-array nvars :element-type 'double-float
                              :initial-contents (sf-l sf)))
         (u (make-array nvars :element-type 'double-float
                              :initial-contents (sf-u sf)))
         (c (make-array nvars :element-type 'double-float
                              :initial-element 0d0))
         (b (coerce (sf-b sf) 'dvec))
         (z0 0d0))
    (map nil (lambda (pair)
               (destructuring-bind (xi . v) pair
                 (setf (aref c xi) v)))
         (sf-c sf))
    (map nil (lambda (triplet)
               (let ((x (triplet-col triplet))
                     (y (triplet-row triplet))
                     (v (triplet-value triplet)))
                 (push (cons x v) (aref constraints y))))
         (sf-A sf))
    (loop for i below ncons
          do (setf (aref constraints i)
                   (make-quadratic (aref constraints i)
                                   :rhs (aref b i))))
    (loop for q across constraints
          for l across lambda
          when (quadratic-p q)
            do (incf z0 (* l (quad-rhs q)))
               (map nil (lambda (i v)
                          (incf (aref c i) (* l v)))
                    (quad-indices q)
                    (quad-coefs q))
               (setf (quad-scale q) (* (sqrt weight) (quad-scale q))))
    (setf (aref constraints ncons)
          (make-linear (loop for i upfrom 0
                             for c across c
                             when (/= c 0d0)
                               collect (cons i c))))
    (let ((constraints (remove nil constraints)))
      (make-approx-state
       (sf-nvars sf) (sf-ncons sf)
       constraints
       nvars
       l u
       weight lambda (- z0)))))

(defun total-violation (subproblem x)
  (let ((acc 0d0))
    (map nil (lambda (c)
               (when (quadratic-p c)
                 (incf acc (expt (violation c x nil) 2))))
         (approx-constraints subproblem))
    (sqrt acc)))

(defstruct (alm-state
            (:conc-name #:alm-))
  standard-instance
  (mu 10d0)
  (omega (/ 10d0))
  (nu (expt (/ 10d0) .1d0))
  (multipliers (error "Missing arg"))
  (multipliers-l (error "Missing arg"))
  (multipliers-u (error "Missing arg")))

(defun clamp (l x u)
  (declare (type dvec l x u))
  (map 'dvec (lambda (l x u)
               (max l (min x u)))
       l x u))

(defun make-alm (sf &optional (mu 10d0) multipliers)
  (let ((low (make-array (sf-ncons sf)
                         :element-type 'double-float
                         :initial-element double-float-negative-infinity))
        (high (make-array (sf-ncons sf)
                          :element-type 'double-float
                          :initial-element double-float-positive-infinity)))
    (loop for i upfrom 0
          for type across (sf-type sf)
          do (case type
               (< (setf (aref low i) 0d0))
               (> (setf (aref high i) 0d0))))
    (make-alm-state :standard-instance sf
                    :mu (float mu 1d0)
                    :omega (/ 1d0 mu)
                    :nu (expt (/ mu) .1d0)
                    :multipliers
                    (or multipliers
                        (make-array (sf-ncons sf)
                                    :element-type 'double-float
                                    :initial-element 0d0))
                    :multipliers-l low
                    :multipliers-u high)))

(defun alm-iteration (state x &optional precision)
  (let ((approx (make-alm-subproblem (alm-standard-instance state)
                                     (alm-multipliers state)
                                     (alm-mu state))))
    (multiple-value-bind (x pg)
        (approx approx 1000000 x
                (or precision
                    (max (alm-omega state) 1d-5)))
      (let ((value (dual-value approx x))
            (violation (map 'dvec (lambda (c)
                                    (if (quadratic-p c)
                                        (violation c x nil)
                                        0d0))
                            (approx-constraints approx))))
        (format t "~12,5g ~12,5g ~g" (2norm violation) pg value)
        (cond ((< (/ (2norm violation) 1)
                  (alm-nu state))
               (setf (alm-multipliers state)
                     (let ((weight (alm-mu state)))
                       (map 'dvec (lambda (l v)
                                    (+ l (* weight v)))
                            (alm-multipliers state) violation))
                     (alm-nu state) (/ (alm-nu state)
                                       (expt (alm-mu state) 0.9d0))
                     (alm-omega state) (max (/ (alm-omega state)
                                               (alm-mu state))
                                            1d-5))
               (format t " minor ~%"))
              (t
               (setf (alm-multipliers state)
                     (let ((weight (alm-mu state)))
                       (map 'dvec (lambda (l v)
                                    (+ l (* weight v)))
                            (alm-multipliers state) violation))
                     (alm-mu state) (min (* 1.5d0 (alm-mu state))
                                         1d6)
                     (alm-nu state) (/ (expt (alm-mu state) 0.1d0))
                     (alm-omega state) (max (/ (alm-mu state))
                                            1d-5))
               (format t " major ~%")))
        (values x violation value)))))

(defun alm-iteration2 (state x &optional precision)
  (let* ((approx (make-alm-subproblem (alm-standard-instance state)
                                      (alm-multipliers state)
                                      (alm-mu state)))
         (violation0 (and x (2norm (map 'dvec (lambda (c)
                                                (if (quadratic-p c)
                                                    (violation c x nil)
                                                    0d0))
                                        (approx-constraints approx))))))
    (multiple-value-bind (x pg)
        (approx approx 1000000 x
                (or precision
                    (max (alm-omega state) 1d-6)))
      (let* ((value (dual-value approx x))
             (violation (map 'dvec (lambda (c)
                                     (if (quadratic-p c)
                                         (violation c x nil)
                                         0d0))
                             (approx-constraints approx)))
             (improvement (and violation0
                               (/ (2norm violation)
                                  violation0))))
        (setf (alm-multipliers state)
              (clamp (alm-multipliers-l state)
                     (let ((weight (alm-mu state)))
                       (map 'dvec (lambda (l v)
                                    (+ l (* weight v)))
                            (alm-multipliers state) violation))
                     (alm-multipliers-u state))
              (alm-mu state) (let ((mu (alm-mu state)))
                               (min (* mu
                                       (if improvement
                                           (max 1d0
                                                (min (* 2d0 improvement)
                                                     2d0))
                                           1d0))
                                    1d7))
              (alm-nu state) (/ (expt (alm-mu state) 0.1d0))
              (alm-omega state) (max (/ (alm-mu state))
                                     1d-6))
        (format t "~12,5g ~12,5g ~12,5g ~g (~12,5g)"
                (inf-norm violation) (2norm violation) pg
                value
                (alm-mu state))
        (values x violation pg value)))))

(defun alm (state &key x0 maxiter)
  (loop with *approx-iterations* = 0
        with x = x0
        with v = nil
        with pg = nil
        with z = nil
        with accuracy = double-float-positive-infinity
        for i below (or maxiter 10000)
        do (setf (values x v pg z) (alm-iteration2
                                    state
                                    x
                                    (min accuracy
                                         (alm-omega state)))
                 v (inf-norm v))
           (setf accuracy (min accuracy (max 1d-5 v)))
           (when (< v 1d-5)
             (setf accuracy 1d-5))
        while (or (> v 1d-5)
                  (and (> pg 1d-5)
                       #+nil(> pg (* 1d-6 (1+ (abs z))))))
        finally (return (values i *approx-iterations* v pg
                                z
                                x))))

(defun next-extrapolation (weight)
  (* .5d0 (1+ (sqrt (1+ (* 4 weight weight))))))

(defun extrapolate (weight prev accelerated current)
  (let* ((next (next-extrapolation weight))
         (vanilla-momentum (/ (1- weight) next))
         (accelerated-momemtun (/ weight next)))
    (map 'dvec (lambda (prev accelerated current)
                 (+ current
                    (* vanilla-momentum
                       (- current prev))
                    (* accelerated-momemtun
                       (- current accelerated))))
         prev accelerated current)))

;; accelerated version... not very good.
(defun aalm (state &key x0 maxiter)
  (loop with *approx-iterations* = 0
        with x = x0
        with v = nil
        with pg = nil
        with z = nil
        with accuracy = double-float-positive-infinity
        with prev-multipliers = (alm-multipliers state)
        for i below (or maxiter 10000)
        for prev-accelerated = (alm-multipliers state)
        for extrapolation = 1d0 then (next-extrapolation extrapolation)
        do (setf (values x v pg z) (alm-iteration2
                                    state
                                    x
                                    (min accuracy
                                         (alm-omega state)))
                 v (inf-norm v))
           (setf accuracy (min accuracy (max 1d-6 v)))
           (when (< v 1d-5)
             (setf accuracy 1d-6))
           (shiftf prev-multipliers
                   (alm-multipliers state)
                   (extrapolate extrapolation
                                prev-multipliers
                                prev-accelerated
                                (alm-multipliers state)))
        while (or (> v 1d-5)
                  (and (> pg 1d-5)
                       (> pg (* 2d-6 (1+ (abs z))))))
        finally (return (values i *approx-iterations* v pg
                                z
                                x))))

(defun adcd-iteration (state x)
  (let* ((approx (make-alm-subproblem (alm-standard-instance state)
                                      (alm-multipliers state)
                                      (alm-mu state)))
         (close (and x
                     (< (2norm
                         (map 'dvec (lambda (c)
                                      (if (quadratic-p c)
                                          (violation c x nil)
                                          0d0))
                              (approx-constraints approx)))
                        5d-2))))
    (multiple-value-bind (x pg)
        (approx approx
                (if close 10000 100)
                x 1d-2)
      (let ((value (dual-value approx x))
            (violation (map 'dvec (lambda (c)
                                    (if (quadratic-p c)
                                        (violation c x nil)
                                        0d0))
                            (approx-constraints approx))))
        (format t "~12,5g ~12,5g ~g" (2norm violation) pg value)
        (let ((out-close (not (not (< pg 5d-2))))
              (almost-approx (not (not (< (2norm violation) 5d-2)))))
          (when (and (< pg 1d-2)
                     (< (2norm violation) 1d-2))
            (throw 'done x))
          (setf (alm-multipliers state)
                (let ((weight (* (if out-close 1d0 .5d0)
                                 (alm-mu state))))
                  (map 'dvec (lambda (l v)
                               (+ l (* weight v)))
                       (alm-multipliers state) violation)))
          (setf (alm-mu state) (min 1d6
                                    (* (cond ((and out-close almost-approx)
                                              1d0)
                                             (out-close
                                              10d0)
                                             (t
                                              1.00d0))
                                       (alm-mu state)))
                (alm-nu state) (/ (expt (alm-mu state) 0.1d0))
                (alm-omega state) (/ (alm-mu state))))
        (values x violation)))))
