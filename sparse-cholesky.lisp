(load-shared-object "/Users/pkhuong/SuiteSparse/libsuitesparse.dylib")

(define-alien-type cholmod-common (struct cholmod-common-struct))

(define-alien-routine cholmod-allocate (* cholmod-common))
(define-alien-routine cholmod-release void (ptr (* cholmod-common)))

(macrolet ((def (suffix type)
             `(progn
                (define-alien-routine ,(alexandria:format-symbol
                                        *package*
                                        "CHOLMOD-GET-~A" suffix)
                  ,type (ptr (* cholmod-common)))
                (define-alien-routine ,(alexandria:format-symbol
                                        *package*
                                        "CHOLMOD-SET-~A" suffix)
                  ,type
                  (ptr (* cholmod-common))
                  (value ,type)))))
  (def #:print int)
  (def #:print-function (* t))
  (def #:dbound double)
  (def #:supernodal-switch double)
  (def #:supernodal int)
  (def #:selected int)
  (def #:itype int) ;; int: 0, long: 2
  (def #:dtype int) ;; double: 0, single: 1
  (def #:status int)
  (def #:fl double)
  (def #:lnz double)
  (def #:anz double)
  (def #:modfl double)
  (def #:malloc-count size-t)
  (def #:memory-usage size-t)
  (def #:memory-inuse size-t)
  (def #:rowfacfl double)
  (def #:aatfl double)
  (def #:blas-ok int))

(define-alien-routine cholmod-start int (ptr (* cholmod-common)))
(define-alien-routine cholmod-finish int (ptr (* cholmod-common)))
(define-alien-routine cholmod-defaults int (ptr (* cholmod-common)))
(define-alien-routine cholmod-free-work int (ptr (* cholmod-common)))

(define-alien-type cholmod-sparse
  (struct cholmod-sparse-struct
          (nrow size-t)
          (ncol size-t)
          (nzmax size-t)

          (p (* t))
          (i (* t))
          (nz (* t))

          (x (* t))
          (z (* t))

          (stype int)
          (itype int)
          (xtype int)
          (dtype int)
          (sorted int)
          (packed int)))

(define-alien-routine cholmod-allocate-sparse (* cholmod-sparse)
  (nrow size-t)
  (ncol size-t)
  (nzmax size-t)
  (sorted int)
  (packed int)
  (stype int)
  (xtype int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-free-sparse int
  (A (* (* cholmod-sparse)))
  (common (* cholmod-common)))

(define-alien-routine cholmod-reallocate-sparse int
  (nznew size-t)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-nnz long
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-transpose cholmod-sparse
  (A (* cholmod-sparse))
  (values int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-transpose-unsym int
  (A (* cholmod-sparse))
  (values int)
  (perm (* int))
  (fset (* int))
  (fsize size-t)
  (F (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-transpose-sym int
  (A (* cholmod-sparse))
  (values int)
  (perm (* int))
  (F (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-ptranspose int
  (A (* cholmod-sparse))
  (values int)
  (perm (* int))
  (fset (* int))
  (fsize size-t)
  (common (* cholmod-common)))

(define-alien-routine cholmod-sort int
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-aat (* cholmod-sparse)
  (A (* cholmod-sparse))
  (fset (* int))
  (fsize size-t)
  (mode int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-copy-sparse (* cholmod-sparse)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-copy (* cholmod-sparse)
  (A (* cholmod-sparse))
  (stype int)
  (mode int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-sparse-xtype int
  (to-xtype int)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-type cholmod-factor (struct cholmod-factor-struct))

(define-alien-routine cholmod-free-factor int
  (L (* (* cholmod-factor)))
  (common (* cholmod-common)))

(define-alien-routine cholmod-pack-factor int
  (L (* cholmod-factor))
  (common (* cholmod-common)))

(define-alien-routine cholmod-copy-factor (* cholmod-factor)
  (L (* cholmod-factor))
  (common (* cholmod-common)))

(define-alien-routine cholmod-factor-xtype int
  (to-xtype int)
  (L (* cholmod-factor))
  (common (* cholmod-common)))

(define-alien-type cholmod-dense
  (struct cholmod-dense-struct
          (nrow size-t)
          (ncol size-t)
          (nzmax size-t)
          (d size-t)
          (x (* t))
          (z (* t))
          (xtype int)
          (dtype int)))

(define-alien-routine cholmod-allocate-dense (* cholmod-dense)
  (nrow size-t)
  (ncol size-t)
  (d size-t)
  (xtype int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-zeros (* cholmod-dense)
  (nrow size-t)
  (ncol size-t)
  (xtype int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-free-dense int
  (X (* (* cholmod-dense)))
  (common (* cholmod-common)))

(define-alien-routine cholmod-sparse-to-dense (* cholmod-dense)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-dense-to-sparse (* cholmod-sparse)
  (X (* cholmod-dense))
  (values int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-copy-dense (* cholmod-dense)
  (X (* cholmod-dense))
  (common (* cholmod-common)))

(define-alien-routine cholmod-dense-xtype int
  (to-xtype int)
  (x (* cholmod-dense))
  (common (* cholmod-common)))

(define-alien-type cholmod-triplet
  (struct cholmod-triplet-struct
          (nrow size-t)
          (ncol size-t)
          (nzmax size-t)
          (nnz size-t)
          (i (* t))
          (j (* t))
          (x (* t))
          (z (* t))
          (stype int)
          (itype int)
          (xtype int)
          (dtype int)))

(define-alien-routine cholmod-allocate-triplet (* cholmod-triplet)
  (nrow size-t)
  (ncol size-t)
  (nzmax size-t)
  (stype int)
  (xtype int)
  (common (* cholmod-common)))

(define-alien-routine cholmod-free-triplet int
  (triplet (* (* cholmod-triplet)))
  (common (* cholmod-common)))

(define-alien-routine cholmod-reallocate-triplet int
  (nznew size-t)
  (triplet (* cholmod-triplet))
  (common (* cholmod-common)))

(define-alien-routine cholmod-sparse-to-triplet (* cholmod-triplet)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-triplet-to-sparse (* cholmod-sparse)
  (triplet (* cholmod-triplet))
  (nzmax size-t)
  (common (* cholmod-common)))

(define-alien-routine cholmod-copy-triplet (* cholmod-triplet)
  (triplet (* cholmod-triplet))
  (common (* cholmod-common)))

(define-alien-routine cholmod-triplet-xtype int
  (to-xtype int)
  (triplet (* cholmod-triplet))
  (common (* cholmod-common)))

(define-alien-routine cholmod-version int
  (version (* int)))

(define-alien-routine cholmod-analyze (* cholmod-factor)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-factorize int
  (A (* cholmod-sparse))
  (L (* cholmod-factor))
  (common (* cholmod-common)))

(define-alien-routine cholmod-solve (* cholmod-dense)
  (sys int)
  (L (* cholmod-factor))
  (B (* cholmod-dense))
  (common (* cholmod-common)))

(define-alien-routine cholmod-solve2 int
  (sys int)
  (L (* cholmod-factor))
  (B (* cholmod-dense))
  (Bset (* cholmod-sparse))

  (X-handle (* (* cholmod-dense)))
  (Xset-handle (* (* cholmod-sparse)))

  (Y-handle (* (* cholmod-dense)))
  (E-handle (* (* cholmod-dense)))

  (common (* cholmod-common)))

(define-alien-routine cholmod-spsolve (* cholmod-sparse)
  (sys int)
  (L (* cholmod-factor))
  (B (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-rcond double
  (L (* cholmod-factor))
  (common (* cholmod-common)))

(define-alien-routine cholmod-print-common int
  (name c-string)
  (common (* cholmod-common)))

(define-alien-routine cholmod-print-sparse int
  (A (* cholmod-sparse))
  (name c-string)
  (common (* cholmod-common)))

(define-alien-routine cholmod-print-dense int
  (A (* cholmod-dense))
  (name c-string)
  (common (* cholmod-common)))

(define-alien-routine cholmod-print-factor int
  (L (* cholmod-factor))
  (name c-string)
  (common (* cholmod-common)))

(define-alien-routine cholmod-print-triplet int
  (truplet (* cholmod-triplet))
  (name c-string)
  (common (* cholmod-common)))

(define-alien-routine cholmod-drop int
  (tol double)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-scale int
  (s (* cholmod-dense))
  (scale int)
  (A (* cholmod-sparse))
  (common (* cholmod-common)))

(define-alien-routine cholmod-sdmult int
  (A (* cholmod-sparse))
  (transpose int)
  (alpha (* double))
  (beta (* double))
  (X (* cholmod-dense))
  (Y (* cholmod-dense))
  (common (* cholmod-common)))

(defvar *cholmod-common*)

(defun make-dense (nrow ncol vector &optional destination)
  (declare (type (simple-array double-float 1) vector))
  (assert (= (length vector) (* nrow ncol)))
  (let ((dense (or destination
                   (cholmod-allocate-dense nrow ncol nrow
                                           1
                                           *cholmod-common*))))
    (declare (type (alien (* cholmod-dense)) dense))
    (assert (eql (slot dense 'dtype) 0))
    (assert (= (slot dense 'nrow) nrow))
    (assert (= (slot dense 'ncol) ncol))
    (sb-kernel:copy-ub8-to-system-area
     vector 0
     (sb-alien:alien-sap (slot dense 'x))
     0
     (* 8 (length vector)))
    dense))

(defun make-dense-from-matlisp (matrix &optional destination)
  (make-dense (matlisp:nrows matrix)
              (matlisp:ncols matrix)
              (matlisp::store matrix)
              destination))

(defun dense-to-matlisp (x &optional destination)
  (declare (type (alien (* cholmod-dense)) x))
  (let* ((nrow (slot x 'nrow))
         (ncol (slot x 'ncol))
         (result (or destination
                     (matlisp:make-real-matrix-dim nrow ncol))))
    (assert (= nrow (matlisp:nrows result)))
    (assert (= ncol (matlisp:ncols result)))
    (sb-kernel:copy-ub8-from-system-area (sb-alien:alien-sap (slot x 'x)) 0
                                         (matlisp::store result)
                                         0
                                         (* 8 nrow ncol))
    result))

(define-alien-routine fflush void (handle (* t)))

(defun flush ()
  (fflush nil))

(defun make-default-common ()
  (let ((common (cholmod-allocate)))
    (cholmod-start common)
    (cholmod-defaults common)
    common))

(defun free-common (common)
  (cholmod-finish common)
  (cholmod-release common)
  nil)

(defmacro with-cholmod (() &body body)
  (let ((temp (gensym "COMMON")))
    `(let* ((,temp (make-default-common))
            (*cholmod-common* ,temp))
       (unwind-protect
            (locally ,@body)
         (free-common ,temp)))))

;; Solve for (A A')x = b
(defun solve-dense (A b &aux (common *cholmod-common*))
  (with-alien ((A (* cholmod-dense) :local (make-dense-from-matlisp A))
               (As (* cholmod-sparse) :local (cholmod-dense-to-sparse
                                              A
                                              1
                                              common))
               (b (* cholmod-dense) :local (make-dense-from-matlisp b))
               (factor (* cholmod-factor) :local (cholmod-analyze As
                                                                  common)))
    (cholmod-set-status common 0)
    (cholmod-factorize As factor common)
    (when (/= (cholmod-get-status common) 0)
      (return-from solve-dense))
    (with-alien ((x (* cholmod-dense) :local (cholmod-solve 0
                                                            factor
                                                            b
                                                            common)))
      (prog1 (dense-to-matlisp x)
        (cholmod-free-dense (addr x) common)
        (cholmod-free-factor (addr factor) common)
        (cholmod-free-dense (addr b) common)
        (cholmod-free-sparse (addr As) common)
        (cholmod-free-dense (addr A) common)))))

(defun make-sparse-from-triplet-vector (nrow ncol vector
                                        &aux (nnz (length vector))
                                          (common *cholmod-common*))
  (with-alien ((triplets (* cholmod-triplet) :local
                         (cholmod-allocate-triplet nrow ncol
                                                   nnz
                                                   0 ;; unsymmetric
                                                   1 ;; real
                                                   common))
               (rows (* int) :local (cast (slot triplets 'i)
                                          (* int)))
               (cols (* int) :local (cast (slot triplets 'j)
                                          (* int)))
               (xs   (* double) :local (cast (slot triplets 'x)
                                             (* double))))
    (loop for i upfrom 0
          for triplet across vector
          do (assert (< (triplet-row triplet) nrow))
             (assert (< (triplet-col triplet) ncol))
             (setf (deref rows i) (triplet-row triplet)
                   (deref cols i) (triplet-col triplet)
                   (deref xs   i) (triplet-value triplet)))
    (setf (slot triplets 'nnz) nnz)
    (let ((sparse (cholmod-triplet-to-sparse triplets nnz common)))
      (cholmod-sort sparse common)
      (cholmod-free-triplet (addr triplets) common)
      sparse)))

(defun scale-sparse! (sparse scale
                     &aux (n (matlisp:nrows scale))
                       (common *cholmod-common*))
  (declare (type (alien (* cholmod-sparse)) sparse))
  (assert (= n (slot sparse 'ncol)))
  (with-alien ((dense (* cholmod-dense) :local
                      (make-dense-from-matlisp scale)))
    (assert (/= (cholmod-scale dense
                               2 ; scale columns
                               sparse common)
                0))
    (cholmod-free-dense (addr dense) common)
    sparse))

(defun scale-sparse (sparse scale)
  (scale-sparse! (cholmod-copy-sparse sparse *cholmod-common*)
                 scale))

(defstruct solve-sparse-state
  factor
  rhs
  solution
  workspace-y
  workspace-e)

(defun free-sparse-state (state)
  (when (solve-sparse-state-factor state)
    (with-alien ((L (* cholmod-factor) :local
                    (shiftf (solve-sparse-state-factor state) nil)))
      (assert (/= 0 (cholmod-free-factor (addr L) *cholmod-common*)))))
  (flet ((free-dense (dense)
           (when dense
             (with-alien ((d (* cholmod-dense) :local
                             dense))
               (assert (/= 0 (cholmod-free-dense (addr d)
                                                 *cholmod-common*)))))))
    (macrolet ((foo (&rest slots)
                 `(progn
                    ,@(loop for slot in slots collect
                            `(free-dense (shiftf (,slot state) nil))))))
      (foo solve-sparse-state-rhs
           solve-sparse-state-solution
           solve-sparse-state-workspace-y
           solve-sparse-state-workspace-e))))

(defun solve-sparse-one-shot (As b &aux (common *cholmod-common*))
  (with-alien ((As (* cholmod-sparse) :local As)
               (b (* cholmod-dense) :local (make-dense-from-matlisp b))
               (factor (* cholmod-factor) :local (cholmod-analyze As
                                                                  common)))
    (cholmod-set-status common 0)
    (cholmod-factorize As factor common)
    (when (/= (cholmod-get-status common) 0)
      (return-from solve-sparse-one-shot))
    (with-alien ((x (* cholmod-dense) :local (cholmod-solve 0
                                                            factor
                                                            b
                                                            common)))
      (prog1 (dense-to-matlisp x)
        (cholmod-free-dense (addr x) common)
        (cholmod-free-factor (addr factor) common)
        (cholmod-free-dense (addr b) common)))))

(defun solve-sparse-recycle (As b state factorized
                             &aux (common *cholmod-common*)
                               (matlisp-b b))
  (with-alien ((As (* cholmod-sparse) :local As)
               (b (* cholmod-dense) :local (make-dense-from-matlisp
                                            b
                                            (solve-sparse-state-rhs state)))
               (factor (* cholmod-factor) :local
                       (or (solve-sparse-state-factor state)
                           (setf (solve-sparse-state-factor state)
                                 (cholmod-analyze As common))))
               (X (* cholmod-dense) :local
                  (solve-sparse-state-solution state))
               (Y (* cholmod-dense) :local
                  (solve-sparse-state-workspace-y state))
               (E (* cholmod-dense) :local
                  (solve-sparse-state-workspace-e state)))
    (unless factorized
      (cholmod-set-status common 0)
      (cholmod-factorize As factor common)
      (when (/= (cholmod-get-status common) 0)
        (return-from solve-sparse-recycle)))
    (unless (plusp (cholmod-solve2
                    0
                    factor
                    b nil
                   
                    (addr X) nil
                    (addr Y) (addr E)

                    common))
      (format t "status: ~A~%" (cholmod-get-status common)))
    (setf (solve-sparse-state-rhs state) b
          (solve-sparse-state-solution state) X
          (solve-sparse-state-workspace-y state) Y
          (solve-sparse-state-workspace-e state) E)
    (dense-to-matlisp x matlisp-b)))

(defun solve-sparse (As b &optional state factorized)
  (if state
      (solve-sparse-recycle As b state factorized)
      (solve-sparse-one-shot As b)))

(defun sparse-m* (sparse x &key transpose
                             y
                             (alpha 1d0)
                             (beta (if y 1d0 0d0))
                             output
                  &aux (nrow (slot sparse 'nrow))
                    (common *cholmod-common*))
  (declare (type (alien (* cholmod-sparse))  sparse)
           (optimize debug))
  (let ((m nrow)
        (n (slot sparse 'ncol)))
    (when transpose (rotatef m n))
    (assert (= 1 (matlisp:ncols x)))
    (assert (= n (matlisp:nrows x)))
    (when y
      (assert (= 1 (matlisp:ncols y)))
      (assert (= m (matlisp:nrows y)))))
  (with-alien ((a (array double 2))
               (b (array double 2))
               (X (* cholmod-dense) :local (make-dense-from-matlisp x))
               (Y (* cholmod-dense) :local
                  (if y
                      (make-dense-from-matlisp y)
                      (cholmod-zeros (if transpose
                                         (slot sparse 'ncol)
                                         nrow)
                                     1
                                     1
                                     common))))
    (setf (deref a 0) alpha
          (deref a 1) 0d0
          (deref b 0) beta
          (deref b 1) 0d0)
    (unless (/= 0 (cholmod-sdmult sparse (if transpose 1 0)
                                  (addr (deref a 0))
                                  (addr (deref b 0))
                                  X
                                  Y
                                  common))
      (flush)
      (cholmod-print-sparse sparse "A" common)
      (cholmod-print-dense X "x" common)
      (cholmod-print-dense Y "y" common)
      (flush)
      (break))
    (prog1 (dense-to-matlisp Y output)
      (cholmod-free-dense (addr Y) common)
      (cholmod-free-dense (addr X) common))))
