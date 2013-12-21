#||

KKT system:

ZSe = 0
WTe = 0
Ax = b
yA + z - w = c

Newton system:

      Z Δx + S Δz      = ZSe
     -W Δx      + T Δw = WTe
      A Δx             = Ax - b
AᵗΔy         + Δz - Δw  = yA + z - w - c   


Block matrix:

Δw      Δx      Δy      Δz    
 U     -W                     = e
        Z               L     = f
        A                     = g
-I              Aᵗ      -I     = h

||#
(defun scale-U (u w e)
  (let ((inv-u (matlisp:map-matrix! #'/ u)))
    (values (matlisp:m.*! inv-u w)
            (matlisp:m.*! inv-u e))))

(defun scale-Z (z l f)
  (scale-U z l f))

#||
After in-place scaling:

Δw      Δx      Δy      Δz    
 I     -W                     = e
        I               L     = f
        A                     = g
-I              Aᵗ       I     = h

Clear Δw column
Δw      Δx      Δy      Δz    
 I     -W                     = e
        I               L     = f
        A                     = g
       -W       Aᵗ       I     = h

Clear Δx column
Δw      Δx      Δy      Δz    
 I     -W                     = e
        I               L     = f
                      -AL     = g
                Aᵗ     WL+I    = h
||#
(defun clear-delta-w (e h)
  (matlisp:m.+! e h))

(defun clear-delta-x-g (A f g)
  (matlisp:gemm! -1 A f 1 g))

(defun clear-delta-x-h (w l f h)
  (values (matlisp:m+! 1 (matlisp:m.* w l))
          (matlisp:m+! (matlisp:m.* w f) h)))
#||
Permutation
Δw      Δx      Δz    Δy
 I     -W                     = e
        I       L             = f
              WL+I    Aᵗ       = h
               -AL            = g

Scale WL+I
Δw      Δx      Δz    Δy
 I     -W                       = e
        I       L               = f
                I     DAᵗ        = h
               -AL              = g

Clear Δz
Δw      Δx      Δz    Δy
 I     -W                       = e
        I       L               = f
                I     DAᵗ        = h
                    ALDAᵗ        = g + ALh = g'

sqrt(LD) = S
||#

(defun scale-wl+1 (wl+1 h)
  (let ((inv (matlisp:map-matrix! #'/ wl+1)))
    (values inv (matlisp:m.*! inv h))))

(defun clear-delta-z (A l d h g)
  (values (matlisp:map-matrix! #'sqrt (matlisp:m.* l d))
          (matlisp:gemm! 1 A (matlisp:m.* l h) 1 g)))

(defvar *worst-condition-number* 1d0)

(defun cond-number (a)
  (multiple-value-bind (eigen ok)
      (matlisp:svd a)
    (when ok
      (let ((cond (/ (matlisp:mref eigen 0)
                     (matlisp:mref eigen (1- (matlisp:nrows eigen))))))
        (setf *worst-condition-number* (max *worst-condition-number*
                                            cond))
        cond))))

(defun solve-delta-y (A S g)
  (let ((diag (matlisp:make-real-matrix-dim (matlisp:nrows S)
                                            (matlisp:nrows S))))
    (setf (matlisp:diag diag) S)
    (solve-dense (matlisp:m* A diag) g)))

#||
Scale ALDAᵗ
Δw      Δx      Δz    Δy
 I     -W                       = e
        I       L               = f
                I     DAᵗ        = h
                      I         = Δy
||#

(defun solve-delta-z (D A h dy)
  (let* ((DAtdy (matlisp:m.*! D
                              (matlisp:m* (matlisp:transpose A)
                                          dy))))
    (matlisp:axpy! -1 DAtdy h)))

(defun solve-delta-x (L f delta-z)
  (matlisp:axpy! -1 (matlisp:m.* delta-z L) f))

(defun solve-delta-w (W e delta-x)
  (matlisp:m+! (matlisp:m.*! delta-x W) e))

(defun solve-kkt-newton (l u w z A
                         e f g h)
  (let (wl+1 d s
        dw dx dy dz)
    (setf (values w e) (scale-U U w e)
          (values l f) (scale-Z Z l f)
          h            (clear-delta-w e h)
          g            (clear-delta-x-g A f g)
          (values wl+1 h) (clear-delta-x-h w l f h)
          (values d h) (scale-wl+1 wl+1 h)
          (values s g) (clear-delta-z A l d h g)
          dy (solve-delta-y A S g)
          dz (solve-delta-z D A H dy)
          dx (solve-delta-x L f dz)
          dw (solve-delta-w W e dx))
    (values dw dx dy dz)))

#||
Δw      Δx      Δy      Δz    
 U     -W                     = e
        Z               L     = f
        A                     = g
-I              Aᵗ       I     = h
||#
(defun delta (expected value)
  (matlisp:norm (matlisp:m- expected value)))

(defun test-kkt-solve (l u w z A e f g h)
  (multiple-value-bind (dw dx dy dz)
      (apply 'solve-kkt-newton
             (mapcar #'matlisp:copy
                     (list l u w z A e f g h)))
    (values (delta (matlisp:m- (matlisp:m.* u dw)
                               (matlisp:m.* w dx))
                   e)
            (delta (matlisp:m+ (matlisp:m.* z dx)
                               (matlisp:m.* l dz))
                   f)
            (delta (matlisp:m* A dx) g)
            (delta (matlisp:m- (matlisp:m+ (matlisp:m* (matlisp:transpose A)
                                                       dy)
                                           dz)
                               dw)
                   h))))

(defun random-positive-vector (n)
  (matlisp:m.+ 1d-1 (matlisp:m.* 10 (matlisp:rand n 1))))

(defun test-m-n (m n)
  (with-cholmod ()
    (test-kkt-solve (random-positive-vector n)
                    (random-positive-vector n)
                    (random-positive-vector n)
                    (random-positive-vector n)

                    (matlisp:m+ (matlisp:rand m n)
                                (matlisp:eye m n))

                    (matlisp:rand n 1)
                    (matlisp:rand n 1)
                    (matlisp:rand m 1)
                    (matlisp:rand n 1))))

(defun test (max &aux (worst 0d0) (*worst-condition-number* 0d0))
  (loop for m from 1 upto max do
    (loop for n from m upto max do
      (loop repeat 10 do
        (let* ((deltas (multiple-value-list (test-m-n m n)))
               (max (reduce #'max deltas)))
          (setf worst (max max worst))
          (when (> max 1d-6)
            (format t "~A ~A ~A~%" m n deltas))))))
  (values worst *worst-condition-number*))
