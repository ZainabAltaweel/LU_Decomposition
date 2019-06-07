namespace LUDecompoistion
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.btnComputeLU = new System.Windows.Forms.Button();
            this.txtMatrixSize = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.txtBlockSize = new System.Windows.Forms.TextBox();
            this.btnLUParallel = new System.Windows.Forms.Button();
            this.SuspendLayout();
            // 
            // btnComputeLU
            // 
            this.btnComputeLU.Location = new System.Drawing.Point(39, 123);
            this.btnComputeLU.Name = "btnComputeLU";
            this.btnComputeLU.Size = new System.Drawing.Size(234, 47);
            this.btnComputeLU.TabIndex = 0;
            this.btnComputeLU.Text = "Compute LU (Sequential)";
            this.btnComputeLU.UseVisualStyleBackColor = true;
            this.btnComputeLU.Click += new System.EventHandler(this.btnComputeLU_Click);
            // 
            // txtMatrixSize
            // 
            this.txtMatrixSize.Location = new System.Drawing.Point(173, 25);
            this.txtMatrixSize.Name = "txtMatrixSize";
            this.txtMatrixSize.Size = new System.Drawing.Size(100, 22);
            this.txtMatrixSize.TabIndex = 1;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(33, 29);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(114, 17);
            this.label1.TabIndex = 2;
            this.label1.Text = "Enter Matrix Size";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(36, 63);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(111, 17);
            this.label2.TabIndex = 3;
            this.label2.Text = "Enter Block Size";
            // 
            // txtBlockSize
            // 
            this.txtBlockSize.Location = new System.Drawing.Point(173, 63);
            this.txtBlockSize.Name = "txtBlockSize";
            this.txtBlockSize.Size = new System.Drawing.Size(100, 22);
            this.txtBlockSize.TabIndex = 4;
            // 
            // btnLUParallel
            // 
            this.btnLUParallel.Location = new System.Drawing.Point(39, 177);
            this.btnLUParallel.Name = "btnLUParallel";
            this.btnLUParallel.Size = new System.Drawing.Size(234, 53);
            this.btnLUParallel.TabIndex = 5;
            this.btnLUParallel.Text = "Compute LU (Parallel)";
            this.btnLUParallel.UseVisualStyleBackColor = true;
            this.btnLUParallel.Click += new System.EventHandler(this.btnLUParallel_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(451, 337);
            this.Controls.Add(this.btnLUParallel);
            this.Controls.Add(this.txtBlockSize);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.txtMatrixSize);
            this.Controls.Add(this.btnComputeLU);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button btnComputeLU;
        private System.Windows.Forms.TextBox txtMatrixSize;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.TextBox txtBlockSize;
        private System.Windows.Forms.Button btnLUParallel;
    }
}

