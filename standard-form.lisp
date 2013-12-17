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
  b
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
                        (vector-push-extend lb b))
                       ((= ub double-float-positive-infinity)
                        ;;      ax >= b
                        ;; -> ax - s = b
                        (vector-push-extend lb b)
                        (artificial-var i -1))
                       ((= lb double-float-negative-infinity)
                        ;;       ax <= b
                        ;; -> ax + s = b
                        (vector-push-extend ub b)
                        (artificial-var i 1))
                       (t
                        ;; l <= ax <= b
                        ;; ->
                        ;; ax + s = b
                        ;;  s <= (b - l)
                        (vector-push-extend ub b)
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
     :l (coerce l 'simple-vector)
     :u (coerce u 'simple-vector)
     :initial-vars nvars)))