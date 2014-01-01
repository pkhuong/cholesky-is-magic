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
      (setf (quad-scale q) (/ norm)))))

(defun violation (quadratic x)
  "Compute ax-b"
  (let ((acc (- (quad-rhs quadratic))))
    (map nil (lambda (i c)
               (incf acc (* c (aref x i))))
         (quad-indices quadratic)
         (quad-coefs quadratic))
    (* acc (quad-scale quadratic))))

;; (x-x0)(y-y0)
;;  or
;; (x0-x)(y-y0)
(defstruct (complementarity
            (:constructor %make-complementarity)
            (:conc-name #:comp-))
  (x-flipped nil :type boolean)
  (x0 0d0 :type double-float)
  (xi (error "Missing arg") :type index)
  (y0 0d0 :type double-float)
  (yi (error "Missing arg") :type index))

;; Global Lipschitz-esque value
(defgeneric accumulate-nu (constraint nu)
  (:method ((c linear) nu)
    ;; Lipschitz constant = 0: we're elementwise
    nu)
  (:method ((c complementarity) nu)
    ;; Lipschitz constant = 0: we're elementwise
    nu)
  (:method ((q quadratic) nu)
    (let ((beta (quad-beta q))
          (scale (quad-scale q)))
      (map nil (lambda (i c)
                 (incf (aref nu i) (* beta
                                      (expt (* c scale) 2))))
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
  (c (error "Missing arg") :type simple-vector)
  (b (error "Missing arg") :type dvec)
  (z0 0d0 :type double-float))

(defun primal-value (approx x)
  (reduce #'+ (approx-c approx)
          :key (lambda (pair)
                 (* (aref x (car pair))
                    (cdr pair)))))

(defun dual-value (approx x)
  (let ((acc 0d0)
        (ncons (approx-orig-cons approx))
        (nvars (approx-orig-vars approx)))
    (loop for b across (approx-b approx)
          for i upfrom nvars
          do (incf acc (* b (aref x i))))
    (loop repeat nvars
          for i upfrom (+ nvars ncons)
          for j upfrom (+ nvars ncons nvars)
          for l across (approx-l approx)
          for u across (approx-u approx)
          for z = (aref x i)
          for w = (aref x j)
          do (when (> z 0)
               (assert (> l double-float-negative-infinity))
               (incf acc (* l z)))
             (when (> w 0)
               (assert (< u double-float-positive-infinity))
               (decf acc (* u w))))
    acc))

(defun complementarity-violation (approx x)
  (let ((acc 0d0)
        (ncons (approx-orig-cons approx))
        (nvars (approx-orig-vars approx)))
    (loop repeat nvars
          for i upfrom (+ nvars ncons)
          for j upfrom (+ nvars ncons nvars)
          for l across (approx-l approx)
          for u across (approx-u approx)
          for z = (aref x i)
          for w = (aref x j)
          for xi across x
          do (when (> z 0)
               (incf acc (* (- xi l) z)))
             (when (> w 0)
               (incf acc (* (- u xi) w))))
    acc))

(defun make-approx-state (orig-vars orig-cons
                          constraints nvars l u c b
                          &key (z0 0d0))
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
                        :c c
                        :b b
                        :z0 z0)))
;; Vars: x, y, z, w
;; quad: |Ax-b|
;; quad: |yA + z - w - c|
;; comp: z(x-l)
;; comp: w(u-x)
(defun make-approx (sf &key complementarity (scale t) (l1-penalty 0d0))
  (let* ((nvars (sf-nvars sf))
         (ncons (sf-ncons sf))
         (n (+ (* 3 nvars) ncons))
         (constraints (make-array (+ 2 n) :initial-element nil))
         (l (make-array n :element-type 'double-float
                          :initial-element double-float-negative-infinity))
         (u (make-array n :element-type 'double-float
                          :initial-element double-float-positive-infinity)))
    (replace l (sf-l sf))
    (replace u (sf-u sf))
    (loop for i below nvars
          for yi upfrom ncons
          for zi upfrom (+ nvars ncons)
          for wi upfrom (+ nvars ncons nvars)
          for li across (sf-l sf)
          for ui across (sf-u sf)
          do (cond ((< li -1d8)
                    (setf (aref l zi) 0d0
                          (aref u zi) 0d0))
                   (t
                    (setf (aref l zi) 0d0)
                    (push (cons zi 1d0) (aref constraints yi))
                    (push (cons zi (- li)) (aref constraints n))
                    (when complementarity
                      (setf (aref constraints zi)
                            (%make-complementarity
                             :x0 li
                             :xi i
                             :yi zi)))))
             (cond ((> ui 1d8)
                    (setf (aref l wi) 0d0
                          (aref u wi) 0d0))
                   (t
                    (setf (aref l wi) 0d0)
                    (push (cons wi -1d0) (aref constraints yi))
                    (push (cons wi ui) (aref constraints n))
                    (when complementarity
                      (setf (aref constraints wi)
                            (%make-complementarity
                             :x-flipped t
                             :x0 ui
                             :xi i
                             :yi wi))))))
    (map nil (lambda (triplet)
               (let ((x (triplet-col triplet))
                     (y (triplet-row triplet))
                     (v (triplet-value triplet)))
                 (push (cons x v) (aref constraints y))
                 (push (cons (+ nvars y) v) (aref constraints
                                                  (+ ncons x)))))
         (sf-A sf))
    (loop for i below ncons
          for yi upfrom nvars
          for pairs across constraints
          for b across (sf-b sf)
          for type across (sf-type sf)
          when pairs
          do (push (cons (+ i nvars) (- b)) (aref constraints n))
             (setf (aref constraints i)
                   (make-quadratic pairs :rhs b))
             (case type
               (< (setf (aref u yi) 0d0))
               (> (setf (aref l yi) 0d0))))
    (loop repeat nvars
          for i from ncons
          do (setf (aref constraints i)
                   (make-quadratic (aref constraints i))))
    (map nil (lambda (pair)
               (destructuring-bind (xi . v) pair
                 (setf (quad-rhs (aref constraints (+ xi ncons))) v)
                 (push (cons xi v) (aref constraints n))))
         (sf-c sf))
    (setf (aref constraints n)
          (make-quadratic (aref constraints n)))
    (setf (aref constraints (1+ n))
          (make-linear
           (loop
             for i below nvars
             for l across l
             for u across u
             collect (cons i
                           (cond ((and (= l double-float-negative-infinity)
                                       (< u double-float-positive-infinity))
                                  (- l1-penalty))
                                 ((and (> l double-float-negative-infinity)
                                       (= u double-float-positive-infinity))
                                  l1-penalty)
                                 (t 0d0)))
             collect (cons (+ i nvars ncons) l1-penalty)
             collect (cons (+ i nvars ncons nvars)
                           l1-penalty))))
    (when scale
      (loop for i upto n
            for c = (aref constraints i)
            when (quadratic-p c)
              do (scale-quadratic c)))
    (let ((constraints (remove nil constraints)))
      (make-approx-state
       (sf-nvars sf) (sf-ncons sf)
       constraints
       n
       l u
       (coerce (sf-c sf) 'simple-vector)
       (coerce (sf-b sf) 'dvec)))))

(defgeneric %value-&-gradient (constraint x &optional g)
  (:method ((c linear) x &optional g)
    (let ((acc 0d0))
      (map nil (lambda (i v)
                 (incf acc (* v (aref x i)))
                 (incf (aref g i) v))
           (linear-indices c)
           (linear-coefs c))
      acc))
  (:method ((c complementarity) x &optional g)
    (let* ((xi (comp-xi c))
           (yi (comp-yi c))
           (xk (- (aref x xi) (comp-x0 c)))
           (yk (- (aref x yi) (comp-y0 c)))
           (flipped (comp-x-flipped c)))
      (when flipped
        (setf xk (- xk)))
      (when (minusp xk)
        (setf xk 0d0))
      (when (minusp yk)
        (setf yk 0d0))
      (when g
        (incf (aref g xi) (if flipped (- yk) yk))
        (incf (aref g yi) xk))
      (* yk xk)))
  (:method ((q quadratic) x &optional g)
    (let ((violation (violation q x))  ;; ax-b
          (scale (quad-scale q)))
      (when g
        (map nil (lambda (i c)
                   (incf (aref g i) (* c scale violation)))
             (quad-indices q)
             (quad-coefs q)))
      (* .5d0 (expt violation 2)))))

(defun value-&-gradient (state x)
  (let ((g (make-array (approx-nvars state)
                       :element-type 'double-float
                       :initial-element 0d0))
        (z 0d0)
        (max 0d0))
    (map nil (lambda (c)
               (let ((v (%value-&-gradient c x g)))
                 (incf z v)
                 (setf max (max max (abs v)))))
         (approx-constraints state))
    (values z
            g
            max)))

(defun solve-coordinate (z nu theta g l u)
  "Min_x gx + (n theta nu)/2tau (x-z)^2, l <= x <= u,
   Tau = n"
  (let ((step (* theta nu)))
    (cond ((zerop step)
           (ecase (signum g)
             (-1d0 u)
             ((0d0 -0d0) z)
             (1d0 l)))
          (t
           (let ((best (- z (/ g step))))
             (cond ((< best l)
                    l)
                   ((> best u)
                    u)
                   (t best)))))))

(defun approx-descent (state theta y z)
  (multiple-value-bind (value g max)
      (value-&-gradient state y)
    (values (map 'dvec
                 (lambda (z nu g l u)
                   (solve-coordinate z nu theta g l u))
                 z
                 (approx-nu state)
                 g
                 (approx-l state)
                 (approx-u state))
            g
            value
            max)))

(defun approx-iteration (state theta x z)
  (let* ((y (map 'dvec (lambda (x z)
                         (+ (* (- 1d0 theta) x)
                            (* theta z)))
                 x z)))
    (multiple-value-bind (zp g)
        (approx-descent state theta y z)
      (values (map 'dvec (lambda (y z zp)
                           (+ y (* theta (- zp z))))
                   y z zp)
              zp
              (* .5d0 (- (sqrt (+ (expt theta 4)
                                  (* 4 (expt theta 2))))
                         (expt theta 2)))
              g))))

(defun project-gradient (state x g)
  (map 'dvec (lambda (l x g u)
               (let ((xp (- x g)))
                 (when (< xp l)
                   (setf xp l))
                 (when (> xp u)
                   (setf xp u))
                 (- x xp)))
       (approx-l state)
       x g
       (approx-u state)))

(defun dot-diff (g x xp)
  (let ((acc 0d0))
    (map nil (lambda (g x xp)
               (incf acc (* g (- xp x))))
         g x xp)
    acc))

(defun project (state x)
  (map 'dvec (lambda (l x u)
               (min u (max l x)))
       (approx-l state)
       x
       (approx-u state)))

(defun approx (state n &optional x)
  (let* ((x (project state
                     (or x
                         (make-array (approx-nvars state)
                                     :element-type 'double-float
                                     :initial-element 0d0))))
         (z (copy-seq x))
         zp gy
         (theta 1d0))
    (dotimes (i n z)
      (setf (values x zp theta gy)
            (approx-iteration state theta x z))
      (multiple-value-bind (value g max)
          (value-&-gradient state zp)
        (cond ((plusp (dot-diff g z zp))
               (format t "Restart ~d~%" (1+ i))
               (setf x z
                  theta 1d0))
              (t
               (setf z zp)))
        (let* ((pg (2norm (project-gradient state z g)))
               (done (< pg 1d-10)))
          (when (or (zerop i)
                    (= i (1- n))
                    (zerop (mod (1+ i) 1000))
                    done)
            (format t "~8d: ~12,5g ~12,5g ~12,5g ~12,5g ~12,5g~%"
                    (1+ i)
                    (2norm g)
                    pg
                    max
                    (+ value (approx-z0 state))
                    (complementarity-violation state z)))
          (when done
            (return z)))))))

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
         (b (coerce (sf-b sf) 'dvec)))
    (map nil (lambda (pair)
               (destructuring-bind (xi . v) pair
                 (setf (aref c xi) v)))
         (sf-c sf))
    (map nil (lambda (triplet)
               (let ((x (triplet-col triplet))
                     (y (triplet-row triplet))
                     (v (triplet-value triplet)))
                 (push (cons x v) (aref constraints y))
                 (incf (aref c x) (* (aref lambda y) v))))
         (sf-A sf))
    (loop for i below ncons
          do (setf (aref constraints i)
                   (make-quadratic (aref constraints i)
                                   :rhs (aref b i))))
    (setf (aref constraints ncons)
          (make-linear (loop for i upfrom 0
                             for c across c
                             when (/= c 0d0)
                               collect (cons i c))))
    (loop for c across constraints
          when (quadratic-p c)
            do (scale-quadratic c)
               (setf (quad-scale c) (* (sqrt weight) (quad-scale c))))
    (let ((constraints (remove nil constraints)))
      (make-approx-state
       (sf-nvars sf) (sf-ncons sf)
       constraints
       nvars
       l u
       (coerce (sf-c sf) 'simple-vector)
       b
       :z0 (- (dot lambda c))))))

