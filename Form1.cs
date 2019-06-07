using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using LUDecomposition;
using System.Diagnostics;
using System.Threading.Tasks;

namespace LUDecompoistion
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void btnComputeLU_Click(object sender, EventArgs e)
        {
            Stopwatch SW = new Stopwatch();
            SW.Start();
            int n = int.Parse(txtMatrixSize.Text);
            int m = int.Parse(txtBlockSize.Text);

            if (n % m > 0)
            {
                //throw new System.ArgumentException("Error\n  Enter a valid block size");
                MessageBox.Show("Error\n  Enter a valid block size");
                return;
            }
            int k = n/m; // block number.
            double[,] L = new double[n, n];
            double[,] U = new double[n, n];
            Matrix m1 = new Matrix(n, n);

            m1[0, 0] = 0.5; m1[0, 1] = 2; m1[0, 2] = 1; m1[0, 3] = 4; m1[0, 4] = 2; m1[0, 5] = 4;
            m1[1, 0] = 1.25; m1[1, 1] = 7; m1[1, 2] = 3.5; m1[1, 3] = 13; m1[1, 4] = 0.5; m1[1, 5] = 1;
            m1[2, 0] = 0.5; m1[2, 1] = 5; m1[2, 2] = 3.5; m1[2, 3] = 10.5; m1[2, 4] = 8; m1[2, 5] = 2;
            m1[3, 0] = 0.5; m1[3, 1] = 6; m1[3, 2] = 6; m1[3, 3] = 19; m1[3, 4] = 3; m1[3, 5] = 6;
            m1[4, 0] = 1; m1[4, 1] = 3; m1[4, 2] = 2; m1[4, 3] = 7; m1[4, 4] = 6; m1[4, 5] = 5;
            m1[5, 0] = 5; m1[5, 1] = 2; m1[5, 2] = 7; m1[5, 3] = 2; m1[5, 4] = 1; m1[5, 5] = 7;
            
            Random rnd = new Random();
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    m1[i, j] = rnd.Next(1, 10);
                }
            }


            // storing the first block into m*m array
            Matrix m11 = new Matrix(m, m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    m11[i, j] = m1[i, j];
            double err = 0;

            double[,] L11I = new double[m, m]; // Lower Inverse of the first block 
            double[,] U11I = new double[m, m]; // Upper Inverse of the first block

            // Step 1 & 2 
            // Decompose A11 into L11 and U11 and store the results in L and U matrix and Compute L11 inverse and U11 inverse in L11I and U11I
            m11.A11(m11, L, U, ref err, L11I, U11I);

            // Step 3 
            // Compute the first row and first column
            m1.Step3(L, U, L11I, U11I);

            // Step 4 
            // Compute the rest of the array 
            m1.Step4(L, U, m);

            SW.Stop();
            MessageBox.Show("Sequential: Time taken = " + SW.ElapsedMilliseconds.ToString());
            // Show results

            Display("Matrix A ", m1, n);
            Display("Lower Matrix ", L, n);
            Display("Upper Matrix ", U, n);

            // verify if LU decomp is correct
            double[,] res = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    for (int x = 0; x < n; x++)
                        res[i, j] += L[i, x] * U[x, j];
            err = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    err += Math.Abs(res[i, j] - m1[i, j]);
            double Er = Math.Round(err, 1, MidpointRounding.AwayFromZero);
            MessageBox.Show("Error = " + Er.ToString());

        }

        private void btnLUParallel_Click(object sender, EventArgs e)
        {
            Stopwatch SW2 = new Stopwatch();
            SW2.Start();
            int n = int.Parse(txtMatrixSize.Text);
            int m = int.Parse(txtBlockSize.Text);

            if (n % m > 0)
            {
                //throw new System.ArgumentException("Error\n  Enter a valid block size");
                MessageBox.Show("Error\n  Enter a valid block size");
                return;
            }
            int k = n / m; // block number.
            double[,] L = new double[n, n];
            double[,] U = new double[n, n];
            Matrix m1 = new Matrix(n, n);

            m1[0, 0] = 0.5; m1[0, 1] = 2; m1[0, 2] = 1; m1[0, 3] = 4; m1[0, 4] = 2; m1[0, 5] = 4;
            m1[1, 0] = 1.25; m1[1, 1] = 7; m1[1, 2] = 3.5; m1[1, 3] = 13; m1[1, 4] = 0.5; m1[1, 5] = 1;
            m1[2, 0] = 0.5; m1[2, 1] = 5; m1[2, 2] = 3.5; m1[2, 3] = 10.5; m1[2, 4] = 8; m1[2, 5] = 2;
            m1[3, 0] = 0.5; m1[3, 1] = 6; m1[3, 2] = 6; m1[3, 3] = 19; m1[3, 4] = 3; m1[3, 5] = 6;
            m1[4, 0] = 1; m1[4, 1] = 3; m1[4, 2] = 2; m1[4, 3] = 7; m1[4, 4] = 6; m1[4, 5] = 5;
            m1[5, 0] = 5; m1[5, 1] = 2; m1[5, 2] = 7; m1[5, 4] = 2; m1[5, 4] = 1; m1[5, 5] = 7;

            Random rnd = new Random();
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    m1[i, j] = rnd.Next(1, 10);
                }
            }
            

            // storing the first block into m*m array
            Matrix m11 = new Matrix(m, m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    m11[i, j] = m1[i, j];
            double err = 0;

            double[,] L11I = new double[m, m]; // Lower Inverse of the first block 
            double[,] U11I = new double[m, m]; // Upper Inverse of the first block

            // Step 1 & 2 
            // Decompose A11 into L11 and U11 and store the results in L and U matrix and Compute L11 inverse and U11 inverse in L11I and U11I
            m11.A11(m11, L, U, ref err, L11I, U11I);

            // Step 3 
            // Compute the first row and first column
            m1.Step3Parallel(L, U, L11I, U11I);

            // Step 4 
            // Compute the rest of the array 
            m1.Step4Parallel(L, U, m);

            SW2.Stop();
            MessageBox.Show("Paralle: Time taken = " + SW2.ElapsedMilliseconds.ToString());
            // Show results

            Display("Matrix A ", m1, n);
            Display("Lower Matrix ", L, n);
            Display("Upper Matrix ", U, n);

            // verify if LU decomp is correct
            double[,] res = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    for (int x = 0; x < n; x++)
                        res[i, j] += L[i, x] * U[x, j];
            err = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    err += Math.Abs(res[i, j] - m1[i, j]);
            double Er = Math.Round(err, 1, MidpointRounding.AwayFromZero);
            MessageBox.Show("Error = " + Er.ToString());
        }
        void Display(string title, Matrix M, int n)
        {
            string out1 = title + "\n";
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    out1 += String.Format("{0:f4}", M[i, j]) + " ";
                }
                out1 += "\n";
            }
            MessageBox.Show(out1);
        }
        void Display(string title, double[,] M, int n)
        {
            string out1 = title + "\n";
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    out1 += String.Format("{0:f4}", M[i, j]) + " ";
                }
                out1 += "\n";
            }
            MessageBox.Show(out1);
        }
    }
}
