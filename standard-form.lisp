;; Standard form:
;;
;;      min c'x
;; s.t.
;;      Ax  = b
;;  l <= x <= u

(defstruct (standard-form
            (:conc-name #:sf-))
  nvars ncons
  c ;; vector of (index . value)
  A ;; vector of triplets
  ;; dense vectors of values
  b type
  l u
  initial-vars)

(defun to-standard-form (mps)
  (declare (type mps-data mps))
  (post-process-mps mps)
  (let* ((triplets (mps-data-triplets mps))
         (ntriplets (length triplets))
         (A (make-array ntriplets
                        :adjustable t
                        :fill-pointer ntriplets
                        :initial-contents triplets))
         (b (make-array (length (mps-data-row-data mps))
                        :element-type 'double-float
                        :adjustable t
                        :fill-pointer 0))
         (type (make-array (length (mps-data-row-data mps))
                           :initial-element nil
                           :adjustable t
                           :fill-pointer 0))
         (nvars (length (mps-data-col-data mps)))
         (l (make-array nvars
                        :element-type 'double-float
                        :adjustable t
                        :fill-pointer 0))
         (u (make-array nvars
                        :element-type 'double-float
                        :adjustable t
                        :fill-pointer 0)))
    (map nil (lambda (col)
               (vector-push-extend (col-data-lb col) l)
               (vector-push-extend (col-data-ub col) u))
         (mps-data-col-data mps))
    (flet ((artificial-var (row coef &optional
                                       (lb 0)
                                       (ub double-float-positive-infinity))
             (assert (<= lb ub))
             (let ((n (length l)))
               (assert (= n (length u)))
               (vector-push-extend (float lb 1d0) l)
               (vector-push-extend (float ub 1d0) u)
               (vector-push-extend (make-triplet :col n
                                                 :row row
                                                 :value (float coef 1d0))
                                   A))))
      (loop for i upfrom 0
            for row across (mps-data-row-data mps)
            do (let ((lb (row-data-lb row))
                     (ub (row-data-ub row)))
                 (cond ((= lb ub)
                        (vector-push-extend nil type)
                        (vector-push-extend lb b))
                       ((= ub double-float-positive-infinity)
                        ;;      ax >= b
                        ;; -> ax - s = b
                        (vector-push-extend lb b)
                        (vector-push-extend '> type)
                        (artificial-var i -1))
                       ((= lb double-float-negative-infinity)
                        ;;       ax <= b
                        ;; -> ax + s = b
                        (vector-push-extend ub b)
                        (vector-push-extend '< type)
                        (artificial-var i 1))
                       (t
                        ;; l <= ax <= b
                        ;; ->
                        ;; ax + s = b
                        ;;  s <= (b - l)
                        (vector-push-extend ub b)
                        (vector-push-extend nil type)
                        (artificial-var i 1 0 (- ub lb)))))))
    (assert (= (length l) (length u)))
    (make-standard-form
     :nvars (length l)
     :ncons (length b)
     :c (let ((c (sort (copy-seq
                        (coerce (mps-data-obj-row mps)
                                'simple-vector))
                       #'< :key #'car)))
          (if (eql (mps-data-sense mps) 'min)
              c
              (map-into c (lambda (x)
                            (cons (car x) (- (cdr x))))
                        c)))
     :A (coerce A 'simple-vector)
     :b (coerce b 'simple-vector)
     :type (coerce type 'simple-vector)
     :l (coerce l 'simple-vector)
     :u (coerce u 'simple-vector)
     :initial-vars nvars)))

(defun rescale-sf (sf)
  (let ((norm #+nil(map '(simple-array double-float 1)
                   (lambda (x)
                     (expt x 2d0))
                   (sf-b sf))
              (make-array (sf-ncons sf)
                          :element-type 'double-float
                          :initial-element 0d0)))
    (map nil (lambda (triplet)
               #+nil(incf (aref norm (triplet-row triplet))
                     (expt (triplet-value triplet) 2))
               (alexandria:maxf (aref norm (triplet-row triplet))
                                (abs (triplet-value triplet))))
         (sf-A sf))
    (map-into norm (lambda (x)
                     (let ((norm x #+nil(sqrt x)))
                       (if (< x 1d-6)
                           1d0
                           (/ norm))))
              norm)
    (let ((b (sf-b sf)))
      (setf (sf-b sf) (map-into b #'* b norm)))
    (map nil (lambda (triplet)
               (setf (triplet-value triplet)
                     (* (triplet-value triplet)
                        (aref norm (triplet-row triplet)))))
         (sf-A sf))
    sf))

#+nil
(defun densify (sf)
  (let ((A (matlisp:make-real-matrix-dim (sf-ncons sf)
                                         (sf-nvars sf))))
    (map nil (lambda (triplet)
               (setf (matlisp:matrix-ref A
                                         (triplet-row triplet)
                                         (triplet-col triplet))
                     (triplet-value triplet)))
         (sf-A sf))
    (make-standard-form
     :nvars (sf-nvars sf)
     :ncons (sf-ncons sf)
     :c (let ((v (matlisp:make-real-matrix-dim (sf-nvars sf) 1)))
          (loop for (i . x) across (sf-c sf)
                do (setf (matlisp:matrix-ref v i) x))
          v)
     :A A
     :b (matlisp:make-real-matrix (sf-b sf))
     :l (matlisp:make-real-matrix (sf-l sf))
     :u (matlisp:make-real-matrix (sf-u sf))
     :initial-vars (sf-initial-vars sf))))
