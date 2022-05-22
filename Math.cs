namespace Linal
{
    class DimensionMismatchException : Exception
    {
        public DimensionMismatchException() { }

        public DimensionMismatchException(string text)
            : base(text) { }
    }

    class Vector
    {
        private readonly double[] _items;

        public int Length { get; }

        public double this[int i]
        {
            get => _items[i];
            set => _items[i] = value;
        }

        public Vector(params double[] list)
        {
            _items = (double[])list.Clone();
            Length = list.GetLength(0);
        }

        public Vector(int length)
        {
            _items = new double[length];
            for (int i = 0; i < length; i++)
                _items[i] = 0;

            Length = length;
        }

        public Vector(Vector v)
        {
            _items = (double[])v._items.Clone();
            Length = v.Length;
        }

        public override string ToString()
        {
            string result = "";

            for (int i = 0; i < Length; i++)
                result += this[i].ToString() + " ";

            return result;
        }

        public static double Dot(Vector a, Vector b)
        {
            if (a.Length != b.Length)
                throw new DimensionMismatchException();

            double result = 0;

            for (int i = 0; i < a.Length; i++)
            {
                result += a[i] * b[i];
            }
            return result;
        }

        #region Operator-overloading helpers

        private bool tryToMultiplyByMatrix(Matrix A)
        {
            if (A.ColsCount != Length)
                return false;

            Vector v = new Vector(this);
            for (int i = 0; i < Length; i++)
            {
                this[i] = Vector.Dot(A.GetRow(i), v);
            }

            return true;
        }

        private bool tryToAdd(Vector v)
        {
            if (v.Length != Length)
                return false;

            for (int i = 0; i < Length; i++)
            {
                _items[i] += v[i];
            }
            return true;
        }

        public override bool Equals(Object? obj)
        {
            if (obj is Vector)
            {
                Vector v = (Vector)obj;
                for (int i = 0; i < v.Length; i++)
                {
                    if (v[i] != this[i])
                        return false;
                }

                return true;
            }

            return false;
        }

        public override int GetHashCode()
        {
            return _items.GetHashCode();
        }

        #endregion

        #region Operator overloadings

        public static Vector operator *(Matrix A, Vector v)
        {
            if (A is null || v is null)
                throw new ArgumentNullException();

            var result = new Vector(v);
            if (result.tryToMultiplyByMatrix(A))
                return result;

            throw new DimensionMismatchException();
        }

        public static Vector operator *(double s, Vector v)
        {
            if (v is null)
                throw new ArgumentNullException();

            var result = new Vector(v);
            for (int i = 0; i < v.Length; i++)
            {
                result[i] *= s;
            }

            return result;
        }

        public static Vector operator +(Vector a, Vector b)
        {
            if (a is null || b is null)
                throw new ArgumentNullException();

            var result = new Vector(a);
            if (result.tryToAdd(b))
                return result;

            throw new DimensionMismatchException();
        }

        public static bool operator ==(Vector v1, Vector v2)
        {
            if (v1 is null)
                return v2 is null;

            return v1.Equals(v2);
        }

        public static bool operator !=(Vector v1, Vector v2)
        {
            return !(v1 == v2);
        }

        #endregion
    }

    class Matrix
    {
        private readonly double[,] _items;

        public int RowsCount { get; }
        public int ColsCount { get; }

        public double this[int row, int col]
        {
            get => _items[row, col];
            set => _items[row, col] = value;
        }

        public Matrix(int dimensions, bool isE)
        {
            if (dimensions < 1)
                throw new DimensionMismatchException("dimensions must be greater than zero");

            RowsCount = dimensions;
            ColsCount = dimensions;
            _items = new double[RowsCount, ColsCount];

            if (isE)
            {
                for (int i = 0; i < RowsCount; i++)
                    _items[i, i] = 1;
            }
        }

        public Matrix(Matrix a)
        {
            RowsCount = a.RowsCount;
            ColsCount = a.ColsCount;

            _items = (double[,])a._items.Clone();
        }

        public Matrix(int rows, int cols)
        {
            if (rows < 1 || cols < 1)
                throw new DimensionMismatchException("dimentions must be greater than zero");

            RowsCount = rows;
            ColsCount = cols;
            _items = new double[RowsCount, ColsCount];
        }

        public Matrix(double[,] a)
        {
            if (a == null)
                throw new ArgumentNullException();

            RowsCount = a.GetLength(0);
            ColsCount = a.GetLength(1);
            _items = (double[,])a.Clone();
        }

        public override string ToString()
        {
            string result = "";
            for (int i = 0; i < RowsCount; i++)
            {
                for (int j = 0; j < ColsCount; j++)
                    result += this[i, j].ToString() + "\t";

                result += "\n";
            }

            return result;
        }

        public bool IsSquare => RowsCount == ColsCount;

        #region Operator-overloading helpers

        private bool tryToAdd(Matrix m)
        {
            if (RowsCount != m.RowsCount || ColsCount != m.ColsCount)
                return false;

            for (int i = 0; i < RowsCount; i++)
                for (int j = 0; j < ColsCount; j++)
                    this[i, j] += m[i, j];

            return true;
        }

        private bool tryToSubtract(Matrix m)
        {
            if (RowsCount != m.RowsCount || ColsCount != m.ColsCount)
                return false;

            for (int i = 0; i < RowsCount; i++)
                for (int j = 0; j < ColsCount; j++)
                    this[i, j] -= m[i, j];

            return true;
        }

        public void Scale(double s)
        {
            for (int i = 0; i < RowsCount; i++)
                for (int j = 0; j < ColsCount; j++)
                    this[i, j] *= s;
        }

        private bool tryToMultiply(Matrix m)
        {
            if (RowsCount != m.ColsCount)
                return false;

            for (int i = 0; i < RowsCount; i++)
                for (int j = 0; j < ColsCount; j++)
                {
                    Vector a = GetRow(i);
                    Vector b = m.GetCol(j);

                    this[i, j] = Vector.Dot(a, b);
                }

            return true;
        }


        #endregion

        public Vector GetRow(int index)
        {
            List<double> items = new List<double>();
            for (int i = 0; i < ColsCount; i++)
            {
                items.Add(this[index, i]);
            }

            return new Vector(items.ToArray());
        }

        public Vector GetCol(int index)
        {
            List<double> items = new List<double>();
            for (int i = 0; i < RowsCount; i++)
            {
                items.Add(_items[i, index]);
            }

            return new Vector(items.ToArray());
        }

        public void SetRow(int index, Vector v)
        {
            if (ColsCount != v.Length)
                throw new DimensionMismatchException();

            for (int i = 0; i < ColsCount; i++)
                this[index, i] = v[i];
        }

        public void SetCol(int index, Vector v)
        {
            if (RowsCount != v.Length)
                throw new DimensionMismatchException();

            for (int i = 0; i < RowsCount; i++)
                this[i, index] = v[i];
        }

        #region Gauss elimination helpers

        private void addRowToRowScaled(int fromIndex, int toIndex, double scalar)
        {
            for (int i = 0; i < ColsCount; i++)
                this[toIndex, i] = scalar * this[fromIndex, i] + this[toIndex, i];
        }

        private void swapRows(int fromIndex, int toIndex)
        {
            Vector v1 = GetRow(fromIndex);
            Vector v2 = GetRow(toIndex);

            SetRow(toIndex, v1);
            SetRow(fromIndex, v2);
        }

        private void multiplyRow(int index, double scalar)
        {
            for (int i = 0; i < ColsCount; i++)
                this[index, i] *= scalar;
        }

        private void divideRow(int index, double scalar)
        {
            for (int i = 0; i < ColsCount; i++)
                this[index, i] /= scalar;
        }

        #endregion

        public void Triangulate()
        {
            for (int i = 0; i < RowsCount; i++)
            {
                for (int j = 0; j < ColsCount; j++)
                {
                    if (j >= i)
                        continue;

                    if (Math.Abs(this[j, j]) < double.Epsilon)
                    {
                        for (int k = i; k < RowsCount; k++)
                        {
                            if (Math.Abs(this[k, i]) > Math.Abs(this[j, j]))
                                swapRows(i, k);
                        }
                    }

                    var x = -(this[i, j] / this[j, j]);

                    if (double.IsNaN(x) || double.IsInfinity(x))
                        SetRow(i, new Vector(ColsCount));
                    else
                        addRowToRowScaled(j, i, x);
                }
            }
        }

        public static Matrix Triangular(Matrix M)
        {
            var result = new Matrix(M);
            result.Triangulate();

            return result;
        }

        public Matrix ToTriangular()
        {
            return Triangular(this);
        }

        public void MakeRowEchelon()
        {
            for (int k = 0; k < RowsCount; k++)
            {
                if (Math.Abs(this[k, k]) < double.Epsilon)
                {
                    for (int i = k; i < RowsCount; i++)
                    {
                        if (Math.Abs(this[i, k]) > Math.Abs(this[k, k]))
                            swapRows(k, i);
                    }
                }

                var pivot = this[k, k];
                divideRow(k, pivot);

                for (int i = 0; i < RowsCount; i++)
                {
                    if (i == k || Math.Abs(this[i, k]) < double.Epsilon)
                        continue;

                    var factor = this[i, k];
                    addRowToRowScaled(k, i, -factor);
                }
            }
        }

        public static Matrix RowEchelon(Matrix M)
        {
            var result = new Matrix(M);
            result.MakeRowEchelon();

            return result;
        }

        public Matrix ToRowEchelon(Matrix M)
        {
            return RowEchelon(this);
        }

        public double Determinant()
        {
            if (!IsSquare)
                throw new DimensionMismatchException("determinant defined only for square matrices");

            if (ColsCount == 3)
                return this[0, 0] * this[1, 1] * this[2, 2]
                    + this[0, 1] * this[1, 2] * this[2, 0]
                    + this[0, 2] * this[1, 0] * this[2, 1]
                    - this[0, 2] * this[1, 1] * this[2, 0]
                    - this[0, 0] * this[1, 2] * this[2, 1]
                    - this[0, 1] * this[1, 0] * this[2, 2];

            if (ColsCount == 2)
                return this[0, 0] * this[1, 1]
                    - this[0, 1] * this[1, 0];

            if (ColsCount == 1)
                return this[0, 0];

            var T = this.ToTriangular();
            double det = 1;

            for (int i = 0; i < RowsCount; i++)
                for (int j = 0; j < ColsCount; j++)
                {
                    if (i == j)
                        det *= T[i, j];
                }

            return det;
        }

        public double Minor(int row, int col)
        {
            if (!IsSquare)
                throw new DimensionMismatchException("minor defined only for square matrices");

            Matrix M = new Matrix(RowsCount - 1, ColsCount - 1);

            for (int i = 0, j = 0; i < RowsCount; i++)
            {
                if (i == row)
                    continue;

                for (int k = 0, u = 0; k < ColsCount; k++)
                {
                    if (k == col)
                        continue;

                    M[j, u] = this[i, k];
                    u++;
                }
                j++;
            }

            return M.Determinant();
        }


        public void Transpose()
        {
            Matrix T = new Matrix(this);

            for (int i = 0; i < ColsCount; i++)
            {
                Vector row = T.GetRow(i);
                this.SetCol(i, row);
            }
        }

        public static Matrix Transposed(Matrix M)
        {
            Matrix T = new Matrix(M);
            T.Transpose();

            return T;
        }

        public Matrix ToTransposed()
        {
            return Transposed(this);
        }

        public void Inverse()
        {
            if (!IsSquare)
                throw new DimensionMismatchException("inverse defined only for square matrices");

            var s = 1 / Determinant();

            Adjugate();
            Scale(s);
        }

        public static Matrix Inversed(Matrix M)
        {
            var result = new Matrix(M);
            result.Inverse();

            return result;
        }

        public Matrix ToInversed()
        {
            return Inversed(this);
        }

        public void Adjugate()
        {
            Matrix M = new Matrix(this);

            for (int i = 0; i < RowsCount; i++)
            {
                for (int j = 0; j < ColsCount; j++)
                {
                    this[i, j] = Math.Pow(-1, i + j) * M.Minor(i, j);
                }
            }

            this.Transpose();
        }

        public static Matrix Adjugated(Matrix M)
        {
            var result = new Matrix(M);
            result.Adjugate();

            return result;
        }

        public Matrix ToAdjugated()
        {
            return Adjugated(this);
        }


        public int Rank()
        {
            int result = 0;

            Matrix M = new Matrix(this);
            M.Triangulate();

            Vector z = new Vector(RowsCount);

            for (int i = 0; i < RowsCount; i++)
            {
                Vector v = M.GetRow(i);

                if (v != z)
                    result++;
            }

            return result;
        }


        public static Matrix Augmented(Matrix M, Vector v)
        {
            if (v.Length != M.RowsCount)
                throw new DimensionMismatchException();

            Matrix result = new Matrix(M.RowsCount, M.ColsCount + 1);
            for (int i = 0; i < M.RowsCount; i++)
                for (int j = 0; j < M.ColsCount; j++)
                    result[i, j] = M[i, j];

            for (int i = 0; i < M.RowsCount; i++)
                result[i, result.ColsCount - 1] = v[i];

            return result;
        }

        #region Operator overloadings

        public static Matrix operator +(Matrix A, Matrix B)
        {
            if (A == null || B == null)
                throw new ArgumentNullException();

            var result = new Matrix(A);
            if (result.tryToAdd(B))
                return result;

            throw new DimensionMismatchException();
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A == null || B == null)
                throw new ArgumentNullException();

            var result = new Matrix(A);
            if (result.tryToSubtract(B))
                return result;

            throw new DimensionMismatchException();
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A == null || B == null)
                throw new ArgumentNullException();

            var result = new Matrix(A);
            if (result.tryToMultiply(B))
                return result;

            throw new DimensionMismatchException();
        }

        public static Matrix operator *(double s, Matrix M)
        {
            if (M == null)
                throw new ArgumentNullException();

            var result = new Matrix(M);
            result.Scale(s);

            return result;
        }

        #endregion
    }

    class Solver
    {
        public static Vector GaussJordan(Matrix M, Vector v)
        {
            var A = Matrix.Augmented(M, v);
            A.MakeRowEchelon();

            return A.GetCol(A.ColsCount - 1);
        }

        public static Vector MatrixSolution(Matrix M, Vector v)
        {
            return M.ToInversed() * v;
        }
    }
}
