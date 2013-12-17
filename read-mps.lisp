;; see http://web.mit.edu/lpsolve/doc/mps-format.htm
(require "split-sequence")
(require "alexandria")

(defstruct row-data
  name type rhs range lb ub)

(defstruct col-data
  ;; default: [0, infty)
  name lb ub)

(defstruct triplet
  col row value)

(defstruct mps-data
  name
  ;; default: min
  (sense nil :type (member nil min max))
  ;; name -> index
  (rows (make-hash-table :test #'equalp) :read-only t)
  ;; index -> (name type rhs)
  (row-data (make-array 16 :adjustable t :fill-pointer 0)
   :read-only t)
  ;; vector of (column . weight) pairs
  (obj-row (make-array 16 :adjustable t :fill-pointer 0)
   :read-only t)
  ;; name -> index
  (columns (make-hash-table :test #'equalp) :read-only t)
  ;; index -> (name lower upper)
  (col-data (make-array 16 :adjustable t :fill-pointer 0)
   :read-only t)
  (triplets (make-array 16 :adjustable t :fill-pointer 0)
   :read-only t))

(defvar *mps-stream*)

(defun tokenize-line (line)
  (values (split-sequence:split-sequence #\Space line
                                         :remove-empty-subseqs t)
          (and (plusp (length line))
               (char/= #\Space (aref line 0)))))

(defun get-line ()
  (let ((line (read-line *mps-stream* nil nil)))
    (if line
        (tokenize-line line)
        (values nil t))))

(defun get-section-line ()
  (multiple-value-bind (tokens sectionp)
      (get-line)
    (assert sectionp)
    (values tokens sectionp)))

(defvar *mps-data*)

(defun read-name (line)
  (assert (equalp "name" (first line)))
  (assert (>= (length line) 2))
  (assert (not (mps-data-name *mps-data*)))
  (setf (mps-data-name *mps-data*) (format nil "~{~A~^ ~}" (rest line)))
  (get-section-line))

(defun read-sense (line)
  (assert (equalp "objsense" (first line)))
  (assert (not (mps-data-sense *mps-data*)))
  (multiple-value-bind (tokens sectionp)
      (get-line)
    (assert (not sectionp))
    (assert (= 1 (length tokens)))
    (setf (mps-data-sense *mps-data*)
          (alexandria:eswitch ((first tokens) :test #'string-equal)
            ("max" 'max)
            ("maximize" 'max)
            ("min" 'min)
            ("minimize" 'min))))
  (get-section-line))

(defun read-rows (line)
  (assert (equalp line '("rows")))
  (let ((rows (mps-data-rows *mps-data*))
        (data (mps-data-row-data *mps-data*))
        (nfree 0))
    (flet ((push-row (type name)
             (assert (not (gethash name rows)))
             (let ((type (alexandria:eswitch
                             (type
                              :test #'string-equal)
                           (#\N nil)
                           (#\E '=)
                           (#\G '>=)
                           (#\L '<=))))
               (cond ((null type)
                      (setf (gethash name rows) (- (incf nfree))))
                     (t
                      (setf (gethash name rows) (length data))
                      (vector-push-extend (make-row-data
                                           :name name
                                           :type type
                                           :rhs  nil)
                                          data))))))
      (loop
       (multiple-value-bind (tokens new-section)
           (get-line)
         (when new-section
           (return (values tokens new-section)))
         (assert (= 2 (length tokens)))
         (push-row (first tokens) (second tokens)))))))

(defun mps-float (string)
  (let ((value (let ((*read-default-float-format* 'double-float))
                 (read-from-string string))))
    (when (integerp value)
      (setf value (float value 1d0)))
    (check-type value double-float)
    value))

(defun read-columns (line)
  (assert (equalp line '("columns")))
  (let ((rows (mps-data-rows *mps-data*))
        (columns (mps-data-columns *mps-data*))
        (data (mps-data-col-data *mps-data*))
        (triplets (mps-data-triplets *mps-data*))
        (obj (mps-data-obj-row *mps-data*)))
    (labels ((ensure-column (name)
               (or (gethash name columns)
                   (prog1 (setf (gethash name columns)
                                (length data))
                     (vector-push-extend (make-col-data :name name)
                                         data))))
             (insert-triplet (column row value)
               (let ((column (ensure-column column))
                     (row    (or (gethash row rows)
                                 (error "Unknown row ~A" row)))
                     (value  (mps-float value)))
                 (check-type column unsigned-byte)
                 (check-type value double-float)
                 (cond ((< row -1))
                       ((= row -1)
                        (vector-push-extend (cons column value)
                                            obj))
                       (t
                        (check-type row unsigned-byte)
                        (vector-push-extend (make-triplet :col column
                                                          :row row
                                                          :value value)
                                            triplets))))))
      (loop (multiple-value-bind (tokens sectionp)
                (get-line)
              (when sectionp
                (return (values tokens sectionp)))
              (assert (member (length tokens) '(3 5)))
              (destructuring-bind (column r1 v1 &optional r2 v2)
                  tokens
                (insert-triplet column r1 v1)
                (when r2
                  (insert-triplet column r2 v2))))))))

(defun read-rhs (line)
  (assert (equalp line '("rhs")))
  (let ((rows (mps-data-rows *mps-data*))
        (data (mps-data-row-data *mps-data*)))
    (flet ((add-rhs (row value)
             (let* ((row (or (gethash row rows)
                             (error "Unknown row ~A" row)))
                    (value (mps-float value))
                    (data (if (minusp row)
                              (return-from add-rhs (Format t "~A~%" row))
                              (aref data row))))
               (assert (not (row-data-rhs data)))
               (setf (row-data-rhs data) value))))
      (let ((rhs-name nil))
        (loop
         (multiple-value-bind (tokens sectionp)
             (get-line)
           (when sectionp
             (return (values tokens sectionp)))
           (assert (member (length tokens) '(3 5)))
           (if rhs-name
               (assert (string-equal rhs-name (first tokens)))
               (setf rhs-name (first tokens)))
           (destructuring-bind (r1 v1 &optional r2 v2)
               (rest tokens)
             (add-rhs r1 v1)
             (when r2
               (add-rhs r2 v2)))))))))

(defun read-ranges (line)
  (assert (equalp line '("ranges")))
  (let ((rows (mps-data-rows *mps-data*))
        (data (mps-data-row-data *mps-data*)))
    (flet ((add-range (row value)
             (let* ((row (or (gethash row rows)
                             (error "Unknown row ~A" row)))
                    (value (mps-float value))
                    (data (if (minusp row)
                              (return-from add-range)
                              (aref data row))))
               (assert (not (row-data-range data)))
               (setf (row-data-range data) value))))
      (let ((range-name nil))
        (loop
         (multiple-value-bind (tokens sectionp)
             (get-line)
           (when sectionp
             (return (values tokens sectionp)))
           (assert (member (length tokens) '(3 5)))
           (if range-name
               (assert (string-equal range-name (first tokens)))
               (setf range-name (first tokens)))
           (destructuring-bind (r1 v1 &optional r2 v2)
               (rest tokens)
             (add-range r1 v1)
             (when r2
               (add-range r2 v2)))))))))

(defun read-bounds (line)
  (assert (equalp line '("bounds")))
  (let ((columns (mps-data-columns *mps-data*))
        (data (mps-data-col-data *mps-data*)))
    (flet ((add-bound (column type value)
             (let* ((column (or (gethash column columns)
                                (error "Unknown column ~A" column)))
                    (value (and value
                                (mps-float value)))
                    (data (if (minusp column)
                              (return-from add-bound)
                              (aref data column))))
               (alexandria:eswitch (type :test #'string-equal)
                 ("LO"
                  (check-type value double-float)
                  ;; (assert (not (col-data-lb data)))
                  (setf (col-data-lb data) value))
                 ("UP"
                  (check-type value double-float)
                  ;; (assert (not (col-data-ub data)))
                  (setf (col-data-ub data) value))
                 ("FX"
                  (check-type value double-float)
                  ;; (assert (not (col-data-lb data)))
                  ;; (assert (not (col-data-ub data)))
                  (setf (col-data-lb data) value
                        (col-data-ub data) value))
                 ("FR"
                  ;; (assert (not (col-data-lb data)))
                  ;; (assert (not (col-data-ub data)))
                  (setf (col-data-lb data) double-float-negative-infinity
                        (col-data-ub data) double-float-positive-infinity))
                 ("MI"
                  ;; (assert (not (col-data-lb data)))
                  ;; (assert (not (col-data-ub data)))
                  (setf (col-data-lb data) double-float-negative-infinity
                        (col-data-ub data) 0d0))
                 ("PL"
                  ;; (assert (not (col-data-lb data)))
                  ;; (assert (not (col-data-ub data)))
                  (setf (col-data-lb data) 0d0
                        (col-data-ub data) double-float-positive-infinity))))))
      (let ((bound-name nil))
        (loop
         (multiple-value-bind (tokens sectionp)
             (get-line)
           (when sectionp
             (return (values tokens sectionp)))
           (destructuring-bind (type bound column &optional value)
               tokens
             (if bound-name
                 (assert (string-equal bound-name bound))
                 (setf bound-name bound))
             (add-bound column type value))))))))

(defun read-mps (stream)
  (let ((*mps-stream* stream)
        (*mps-data* (make-mps-data)))
    (multiple-value-bind (header sectionp)
        (get-line)
      (assert sectionp)
      (setf header (read-name header))
      (when (equalp header '("objsense"))
        (setf header (read-sense header)))
      (setf header (read-rows header)
            header (read-columns header)
            header (read-rhs header))
      (when (equalp header '("ranges"))
        (setf header (read-ranges header)))
      (when (equalp header '("bounds"))
        (setf header (read-bounds header)))
      (assert (equalp header '("endata"))))
    *mps-data*))

(defun read-mps-file (file)
  (with-open-file (s file :external-format :utf-8)
    (read-mps s)))

(defun post-process-mps (mps)
  (declare (type mps-data mps))
  (unless (mps-data-sense mps)
    (setf (mps-data-sense mps) 'min))
  (map nil (lambda (row)
             (let* ((rhs (or (row-data-rhs row)
                             0d0))
                    (range (row-data-range row))
                    (abs-range (and range (abs range))))
               (check-type rhs double-float)
               (check-type range (or null double-float))
               (setf (values (row-data-lb row)
                             (row-data-ub row))
                     (if range
                         (ecase (row-data-type row)
                           (<= (values (- rhs abs-range) rhs))
                           (>= (values rhs (+ rhs abs-range)))
                           (=  (if (minusp range)
                                   (values (+ rhs range) rhs)
                                   (values rhs (+ rhs range)))))
                         (ecase (row-data-type row)
                           (<= (values double-float-negative-infinity rhs))
                           (>= (values rhs double-float-positive-infinity))
                           (=  (values rhs rhs)))))
               (assert (<= (row-data-lb row) (row-data-ub row)))))
       (mps-data-row-data mps))
  (map nil (lambda (col)
             (unless (col-data-lb col)
               (setf (col-data-lb col) 0d0))
             (unless (col-data-ub col)
               (setf (col-data-ub col) double-float-positive-infinity)))
       (mps-data-col-data mps)))
