using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Diagnostics;

namespace LUDecomposition
{
    class Matrix  // uses inner class to delegate the work to a function
    {             // for launching it as a separate thread and to
                  // pass data to it and collect results from it
        int m;
        public int M
        {
            get { return m; }
        }
        int n;
        public int N
        {
            get { return n; }
        }

        double[,] data = null; //indexer
        public double this[int i, int j]
        {
            get { return data[i,j]; }
            set { data[i,j] = value; }
        }


        public Matrix(int m1, int n1)
        {
            data = new double[m1, n1];
            m = m1;
            n = n1;
        }

        public Matrix(double[,] dt)
        {
            data = dt;
            m = dt.GetLength(0);
            n = dt.GetLength(1);
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            Matrix mres = new Matrix(A.m, B.n);
            Thread[] thArr = new Thread[A.m];
            for (int i = 0; i < A.m; i++)
            {
                //for (int j = 0; j < B.n; j++)
                //    for (int k = 0; k < A.n; k++)
                //        mres.data[i, j] += A.data[i, k] * B.data[k, j];
                MatrixMul mm = new MatrixMul();
                mm.M1 = A;
                mm.M2 = B;
                mm.Mress = mres;
                mm.Iter = i;
                thArr[i] = new Thread(new ThreadStart(mm.ComputeRow));
                thArr[i].Start();

            }

            for (int p = 0; p < A.m; p++)
            {
                thArr[p].Join();  // wait for all threads to finish
            }

            return mres;

        }
        public static Matrix operator +(Matrix A, Matrix B)
        {
            Matrix mres = new Matrix(A.m, B.n);
            for (int i = 0; i < A.m; i++)
                for (int j = 0; j < A.m; j++)
                    mres[i, j] = A[i, j] + B[i, j];
            return mres;
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
            Matrix mres = new Matrix(A.m, B.n);
            for (int i = 0; i < A.m; i++)
                for (int j = 0; j < A.m; j++)
                    mres[i, j] = A[i, j] - B[i, j];
            return mres;
        }
        public void LUDecompose(double[,]L, double[,] U, ref double error)  // operates on M
        {
            if (m != n)
                throw new Exception("width and column dimensions are not the same..");

            // copy data into A matrix
            double[,] A = new double[n, n];
            A = (double [,])data.Clone();

            for (int k = 0; k < n; k++)
            {
                U[k, k] = data[k, k];
                for (int j = k+1; j < n; j++)
                {
                     U[k, j] = data[k, j];
                }
                for (int i = k; i < n; i++)
                {
                    if (i == k)
                        L[i, k] = 1;
                    else
                        L[i, k] = data[i, k] / U[k, k];
                }
                for (int i = k+1; i < n; i++)
                {
                    for (int j = k+1; j < n; j++)
                        data[i, j] = data[i, j] - L[i, k] * U[k, j];
                   
                }
                
            }

            // verify if LU decomp is correct
            double [,] res = new double[n,n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    for (int k = 0; k < n; k++)
                        res[i, j] += L[i, k] * U[k, j];
            error = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    error += Math.Abs(res[i, j] - A[i, j]);
        }
        public double[,] LInverse()  // assumes lower triangular matrix
        {   // with 1 in diagonals
            // see handout on computing inverses of triangular matrices on CS590 web site
            // the following code implements the algorithm described in the handout
            if (m != n)
                throw new Exception("Matrix is not square..");
            double[,] L = (double[,])data.Clone();
            for (int i = 0; i < n; i++)
            {
                for (int j = i - 1; j >= 0; j--)
                {
                    double sum = 0;
                    if (i == (j + 1))  // first diagonal below main diagonal
                    {
                        for (int k = 0; k < i; k++)
                        {
                            sum += L[i, k] * L[k, j];
                        }
                    }
                    else if (i > (j+2))  //  third diagonal and below
                    {
                        for (int k = j + 1; k <= i; k++)
                        {
                            sum += L[i, k] * data[k, j];
                        }
                    }
                    else if (i > (j+1))  // second diagonal
                    {
                        for (int k = j + 1; k < i; k++)
                        {
                            sum += L[i, k] * L[k, j];
                        }
                    }
                    if (i == (j + 1))
                        L[i, j] = 0 - sum;
                    else if (i > (j +2))
                        L[i,j]= -1 * sum;
                    else if (i > (j + 1))
                        L[i, j] = sum - data[i, j];
                    else
                        L[i, j] = -1*sum ;
                }

            }
            return L;  // L contains inverse of L
        }

        
        public double[,] UInverse()  // assumes upper triangular matrix
        {  // see handout on computing inverses of triangular matrices on CS590 web site
            if (m != n)
                throw new Exception("Matrix is not square..");
            double[,] U = (double[,])data.Clone();
            // find the transpose of the matrix
            double[,] UT = this.Transpose();   // convert U to L
            double[,] D = new double[m, n];  // diagonal matrix
            for (int i = 0; i < m; i++)      
                D[i, i] = 1 / UT[i, i];  // D contains Dinverse
            
            double [,] C = new double[m,n];  // compute C = Dinverse * transpose of U
            for (int i = 0; i < n; i++)      // to make diagonal entries of C to be 1
                for (int j = 0; j < n; j++)
                    for (int k = 0; k < n; k++)
                        C[i, j] +=  D[i,k]* UT[k, j] ;  // now C is like transpose of U but diagonal 1's

            double[,] L = new double[n, n];   // copy C into L. We already know how to do Linverse
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    L[i, j] = C[i, j];
            //----compute Linv, then multiply by Dinv
            for (int i = 0; i < n; i++)
            {
                for (int j = i - 1; j >= 0; j--)
                {
                    double sum = 0;
                    if (i == (j + 1))  // first diagonal below main diagonal
                    {
                        for (int k = 0; k < i; k++)
                        {
                            sum += L[i, k] * L[k, j];
                        }
                    }
                    else if (i > (j + 2))  //  third diagonal and below
                    {
                        for (int k = j + 1; k <= i; k++)
                        {
                            sum += L[i, k] * C[k, j];
                        }
                    }
                    else if (i > (j + 1))  // second diagonal
                    {
                        for (int k = j + 1; k < i; k++)
                        {
                            sum += L[i, k] * L[k, j];
                        }
                    }
                    if (i == (j + 1))
                        L[i, j] = 0 - sum;
                    else if (i > (j + 2))
                        L[i, j] = -1 * sum;
                    else if (i > (j + 1))
                        L[i, j] = sum - C[i, j];
                    else
                        L[i, j] = -1 * sum;
                }

            }
 
            double[,] Res = new double[n, n];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                        Res[i,j] += L[i, k] * D[k,j];  // multiply L which is Cinverse by Dinverse 
                }
                
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    U[i, j] = Res[j, i];  // transpose to obtain Uinverse
                }
            }
            return U;  // U now contains inverse of U
        }

        public double[,] Transpose()
        {
            // transposes the matrix i.e, rows become columns and columns become rows
            double[,] T = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    T[i, j] = data[j, i];
                }
            }
            return T;
        }


        public override string ToString()
        {
            string out1 = "";
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                    out1 += data[i, j].ToString() + " \t";
                out1 += "\n";
            }
            return out1;
        }

        public void A11(Matrix m11, double[,] L, double[,] U, ref double error, double[,] LI, double[,] UI)
        {
            int BlockSize = m11.m;
            // Decompose m11
            m11.LUDecompose(L, U, ref error);

            // Compute L11 inverse and U11 inverse 
            Matrix L11 = new Matrix(L);
            Matrix U11 = new Matrix(U);
            double[,] L11I = new double[BlockSize, BlockSize];
            double[,] U11I = new double[BlockSize, BlockSize];
            L11I = L11.LInverse();
            U11I = U11.UInverse();

            for (int i = 0; i < BlockSize; i++)
                for (int j = 0; j < BlockSize; j++)
                {
                    LI[i, j] = L11I[i, j];
                    UI[i, j] = U11I[i, j];
                }
        }

        public void Step3(double[,] L, double[,] U, double[,] L11I, double[,] U11I)
        {
            int n = L11I.GetLength(0); // block size 
            int k = m / n; //No of blocks

            Matrix L1I = new Matrix(L11I);
            Matrix U1I = new Matrix(U11I);

            Matrix LowerBlock = new Matrix(n, n);
            Matrix UpperBlock = new Matrix(n, n);
            Matrix LowerMul = new Matrix(n, n);
            Matrix UpperMul = new Matrix(n, n);
            object olock = new object();

            //Parallel.For(1, k, p =>
            for (int p = 1; p < k; p++)
            {
                for (int i = 0; i < n; i++)     // Copying the content of each block and store it in a matrix 
                    for (int j = 0; j < n; j++)
                    {
                        LowerBlock[i, j] = data[p * n + i, j];
                        UpperBlock[i, j] = data[i, p * n + j];
                    }
                for (int i = 0; i < n; i++)      // copy the column to L
                    for (int j = 0; j < n; j++)
                    {
                        L[p * n + i, j] = (LowerBlock * U1I)[i, j];
                        U[i, p * n + j] = (L1I * UpperBlock)[i, j];
                    }
            }
        }

        public void Step3Parallel(double[,] L, double[,] U, double[,] L11I, double[,] U11I)
        {
            int n = L11I.GetLength(0); // block size 
            List<Matrix> UpperList = new List<Matrix>();
            List<Matrix> LowerList = new List<Matrix>();

            int k = m / n; //No of blocks
            object olock = new object();
            for (int p = 1; p < k; p++)
            {
                Matrix LowerBlock = new Matrix(n, n);
                Matrix UpperBlock = new Matrix(n, n);
                for (int i = 0; i < n; i++)     // Copying the content of each block and store it in a matrix 
                    for (int j = 0; j < n; j++)
                    {
                        LowerBlock[i, j] = data[p * n + i, j];
                        UpperBlock[i, j] = data[i, p * n + j];
                    }
                LowerList.Add(LowerBlock);
                UpperList.Add(UpperBlock);
            }
            Matrix L1I = new Matrix(L11I);
            Matrix U1I = new Matrix(U11I);
            Parallel.For(1, k, (p) =>
            //for (int p = 1; p < k; p++)
            {
                Parallel.Invoke(
                    () => Column(L, LowerList.ElementAt(p - 1), U1I, p, n),
                    () => Row(U, L1I, UpperList.ElementAt(p - 1), p, n)
                    );
            });
        }

        void Column(double[,] L, Matrix B1, Matrix B2, int p, int n)
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    L[p * n + i, j] = (B1 * B2)[i, j];
                }
        }
        void Row(double[,] U, Matrix B1, Matrix B2, int p, int n)
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    U[i, p * n + j] = (B1 * B2)[i, j];
                }
        }

        public void Step4(double[,] L, double[,] U, int BlockSize)
        {
            n = BlockSize;

            Matrix LBlock = new Matrix(n, n);
            Matrix UBlock = new Matrix(n, n);
            Matrix Mul = new Matrix(n, n);
            Matrix Block = new Matrix(n, n);
            Matrix Prime = new Matrix(n, n);
            double[,] Lower = new double[n, n];
            double[,] Upper = new double[n, n];
            Matrix Sum = new Matrix(n, n);
            double err = 0;
            object olock = new object();
            int k = m / n; //No of blocks
            for (int q = 1; q < k; q++)
            {
                // Process the diagonal 
                for (int s = 0; s < q; s++)
                {
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++)
                        {
                            LBlock[i, j] = L[q * n + i, (s * n) + j];
                            UBlock[i, j] = U[(s * n) + i, q * n + j];
                            Block[i, j] = data[q * n + i, q * n + j];

                        }
                    //MessageBox.Show("p= " + q + "  s=  " + s + "  Uper Block \n" + UBlock );
                    //MessageBox.Show("p= " + q + "  s=  " + s + "  Lower Block \n" + LBlock);

                    Mul = LBlock * UBlock;
                    Sum = Sum + Mul;
                    //MessageBox.Show("p= " + q + "s=  " + s + "Mul \n" + Mul);
                    //MessageBox.Show("p= " + q + "s=  " + s + "sum \n" + Sum);
                }
                Prime = Block - Sum;
                Prime.LUDecompose(Lower, Upper, ref err);

                Matrix UpperBlock = new Matrix(n, n);
                double[,] UpperBlockInverse = new double[n, n];


                Matrix LowerBLock = new Matrix(n, n);
                double[,] LowerBlockInverse = new double[n, n];


                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        UpperBlock[i, j] = Upper[i, j];
                        LowerBLock[i, j] = Lower[i, j];
                    }
                UpperBlockInverse = UpperBlock.UInverse();
                LowerBlockInverse = LowerBLock.LInverse();

                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        UpperBlock[i, j] = UpperBlockInverse[i, j];
                        LowerBLock[i, j] = LowerBlockInverse[i, j];
                    }

                //MessageBox.Show("q= " + q + "  Upper Block digonal\n" + UpperBlock);
                //MessageBox.Show("q= " + q + "  Lower diagonal \n" + LowerBLock);

                for (int p = q+1; p < k; p++)
                {
                    //variable to process the upper diagonal matrix 
                    Matrix L1 = new Matrix(n, n);
                    Matrix U1 = new Matrix(n, n);
                    Matrix B1 = new Matrix(n, n);
                    Matrix Sum1 = new Matrix(n, n);
                    Matrix Mul1 = new Matrix(n, n);
                    Matrix Prime1 = new Matrix(n, n);
                    Matrix Final1 = new Matrix(n, n);

                    // vaible to process the lower diagonal matrix 
                    Matrix L2 = new Matrix(n, n);
                    Matrix U2 = new Matrix(n, n);
                    Matrix B2 = new Matrix(n, n);
                    Matrix Sum2 = new Matrix(n, n);
                    Matrix Mul2 = new Matrix(n, n);
                    Matrix Prime2 = new Matrix(n, n);
                    Matrix Final2 = new Matrix(n, n);

                    int r = Math.Min(p, q);
                    for (int s = 0; s < r; s++)
                    {
                        for (int i = 0; i < n; i++)
                            for (int j = 0; j < n; j++)
                            {
                                L1[i, j] = L[q * n + i, s * n + j];
                                U1[i, j] = U[s * n + i, p * n + j];
                                B1[i, j] = data[q * n + i, p * n + j];

                                L2[i, j] = L[p * n + i, s * n + j];
                                U2[i, j] = U[s * n + i, q * n + j];
                                B2[i, j] = data[p * n + i, q * n + j];
                            }
                        // MessageBox.Show("q= " + q + " s= " + s + "  Upper diagonal proccessing L1 \n" + L1);
                        // MessageBox.Show("q= " + q + " s= " + s + "  Upper diagonal proccessing U1 \n" + U1);
                        // MessageBox.Show("q= " + q + " s= " + s + "  Upper diagonal proccessing B1 \n" + B1);


                        Mul1 = L1 * U1;
                        Mul2 = L2 * U2;
                        Sum1 = Sum1 + Mul1;
                        Sum2 = Sum2 + Mul2;
                        //MessageBox.Show("q= " + q + " p =" + p + " s= " + s + " Sum of upper diagonal \n" + Sum1);
                        //MessageBox.Show("q= " + q + " p =" + p + " s= " + s + " Sum of lower diagonal \n" + Sum2);
                    }
                    Prime1 = B1 - Sum1;
                    Prime2 = B2 - Sum2;

                    Final1 = LowerBLock * Prime1;// upper diagonal
                    Final2 = Prime2 * UpperBlock;//  lower diagonal

                    // MessageBox.Show("q= " + q + " Final Upper Block F1 \n" + Final1);

                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++)
                        {
                            L[p * n + i, q * n + j] = Final2[i, j];
                            U[q * n + i, p * n + j] = Final1[i, j];
                        }
                }
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        U[q * n + i, q * n + j] = Upper[i, j];
                        L[q * n + i, q * n + j] = Lower[i, j];
                        Sum[i, j] = 0;
                    }
            }
        }
        public void Step4Parallel(double[,] L, double[,] U, int BlockSize)
        {
            n = BlockSize;

            Matrix LBlock = new Matrix(n, n);
            Matrix UBlock = new Matrix(n, n);
            Matrix Mul = new Matrix(n, n);
            Matrix Block = new Matrix(n, n);
            Matrix Prime = new Matrix(n, n);
            double[,] Lower = new double[n, n];
            double[,] Upper = new double[n, n];
            Matrix Sum = new Matrix(n, n);
            double err = 0;
            object olock = new object();
            int k = m / n; //No of blocks
            for (int q = 1; q < k; q++)
            {
                // Process the diagonal 
                for (int s = 0; s < q; s++)
                {
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++)
                        {
                            LBlock[i, j] = L[q * n + i, (s * n) + j];
                            UBlock[i, j] = U[(s * n) + i, q * n + j];
                            Block[i, j] = data[q * n + i, q * n + j];

                        }
                    //MessageBox.Show("p= " + q + "  s=  " + s + "  Uper Block \n" + UBlock );
                    //MessageBox.Show("p= " + q + "  s=  " + s + "  Lower Block \n" + LBlock);

                    Mul = LBlock * UBlock;
                    Sum = Sum + Mul;
                    //MessageBox.Show("p= " + q + "s=  " + s + "Mul \n" + Mul);
                    //MessageBox.Show("p= " + q + "s=  " + s + "sum \n" + Sum);
                }
                Prime = Block - Sum;
                Prime.LUDecompose(Lower, Upper, ref err);

                Matrix UpperBlock = new Matrix(n, n);
                double[,] UpperBlockInverse = new double[n, n];


                Matrix LowerBLock = new Matrix(n, n);
                double[,] LowerBlockInverse = new double[n, n];


                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        UpperBlock[i, j] = Upper[i, j];
                        LowerBLock[i, j] = Lower[i, j];
                    }
                UpperBlockInverse = UpperBlock.UInverse();
                LowerBlockInverse = LowerBLock.LInverse();

                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        UpperBlock[i, j] = UpperBlockInverse[i, j];
                        LowerBLock[i, j] = LowerBlockInverse[i, j];
                    }

                //MessageBox.Show("q= " + q + "  Upper Block digonal\n" + UpperBlock);
                //MessageBox.Show("q= " + q + "  Lower diagonal \n" + LowerBLock);

                Parallel.For(q + 1, k, p =>
                //for (int p = q+1; p < k; p++)
                {
                    //variable to process the upper diagonal matrix 
                    Matrix L1 = new Matrix(n, n);
                    Matrix U1 = new Matrix(n, n);
                    Matrix B1 = new Matrix(n, n);
                    Matrix Sum1 = new Matrix(n, n);
                    Matrix Mul1 = new Matrix(n, n);
                    Matrix Prime1 = new Matrix(n, n);
                    Matrix Final1 = new Matrix(n, n);

                    // vaible to process the lower diagonal matrix 
                    Matrix L2 = new Matrix(n, n);
                    Matrix U2 = new Matrix(n, n);
                    Matrix B2 = new Matrix(n, n);
                    Matrix Sum2 = new Matrix(n, n);
                    Matrix Mul2 = new Matrix(n, n);
                    Matrix Prime2 = new Matrix(n, n);
                    Matrix Final2 = new Matrix(n, n);

                    int r = Math.Min(p, q);
                    for (int s = 0; s < r; s++)
                    {
                        for (int i = 0; i < n; i++)
                            for (int j = 0; j < n; j++)
                            {
                                L1[i, j] = L[q * n + i, s * n + j];
                                U1[i, j] = U[s * n + i, p * n + j];
                                B1[i, j] = data[q * n + i, p * n + j];

                                L2[i, j] = L[p * n + i, s * n + j];
                                U2[i, j] = U[s * n + i, q * n + j];
                                B2[i, j] = data[p * n + i, q * n + j];
                            }
                        // MessageBox.Show("q= " + q + " s= " + s + "  Upper diagonal proccessing L1 \n" + L1);
                        // MessageBox.Show("q= " + q + " s= " + s + "  Upper diagonal proccessing U1 \n" + U1);
                        // MessageBox.Show("q= " + q + " s= " + s + "  Upper diagonal proccessing B1 \n" + B1);


                        Mul1 = L1 * U1;
                        Mul2 = L2 * U2;
                        Sum1 = Sum1 + Mul1;
                        Sum2 = Sum2 + Mul2;
                        //MessageBox.Show("q= " + q + " p =" + p + " s= " + s + " Sum of upper diagonal \n" + Sum1);
                        //MessageBox.Show("q= " + q + " p =" + p + " s= " + s + " Sum of lower diagonal \n" + Sum2);
                    }
                    Prime1 = B1 - Sum1;
                    Prime2 = B2 - Sum2;

                    Final1 = LowerBLock * Prime1;// upper diagonal
                    Final2 = Prime2 * UpperBlock;//  lower diagonal

                    // MessageBox.Show("q= " + q + " Final Upper Block F1 \n" + Final1);

                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++)
                        {
                            L[p * n + i, q * n + j] = Final2[i, j];
                            U[q * n + i, p * n + j] = Final1[i, j];
                        }
                });
                for (int i=0; i< n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        U[q * n + i, q * n + j] = Upper[i, j];
                        L[q * n + i, q * n + j] = Lower[i, j];
                        Sum[i, j] = 0;
                    }
            }
        }
        class MatrixMul
        {
            Matrix m1;
            public Matrix M1
            {
                get { return m1; }
                set { m1 = value; }
            }
            Matrix m2;
            public Matrix M2
            {
                get { return m2; }
                set { m2 = value; }
            }

            Matrix mress;
            public Matrix Mress
            {
                get { return mress; }
                set { mress = value; }
            }

            int iter;
            public int Iter
            {
                get { return iter; }
                set { iter = value; }
            }

            public void ComputeRow()
            {
                for (int k = 0; k < M1.n; k++)  // changing order still produces
                for (int j = 0; j < M2.n; j++)  // correct result
                   // for (int k = 0; k < M1.n; k++)
                        mress.data[iter, j] += M1.data[iter, k] * M2.data[k, j];

            }
        }
    }
}
