using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InoueLab;
using CsvFileIO;
using System.Numerics;

namespace GaussianProcessRegression
{
    class Program
    {
        static void Main(string[] args)
        {
            double[][] dataX_Y = arrayTojag(CsvFileIO.CsvFileIO.ReadData("GaussianProcessData20121209.txt"));
            double[][] data_X = new double[dataX_Y.GetLength(0)][];
            for (int i = 0; i < dataX_Y.GetLength(0); i++)
            {
                data_X[i] = new double[1];
            }

            double[] data_y = new double[dataX_Y.GetLength(0)];
            double[][] new_data_X = new double[100][];
            for (int i = 0; i < new_data_X.GetLength(0); i++)
            {
                new_data_X[i] = new double[1];
            }
            double[] exp;  //本当はデータ数×d次元配列
            double[] var;
            double alpha = 0.16;
            double beta = 1.0;
            for (int i = 0; i < data_X.Length; i++)
            {
                data_X[i][0] = dataX_Y[i][0];
                data_y[i] = dataX_Y[i][1];
            }
            for (int i = 0; i < new_data_X.Length; i++)
            {
                new_data_X[i][0] = -7.0 + i * 0.1;
            }

            EstimationDistribution(data_X, data_y, new_data_X, alpha, beta, out exp, out var);

        }

        static void EstimationDistribution(double[][] data_X, double[] data_y, double[][] new_data_X, double alpha, double beta, out double[] exp, out double[] var)
        {
            int I = data_X.Length; //データ数
            int J = data_X[0].Length; //次元数;
            exp = new double[new_data_X.Length]; //本当はデータ数×d次元配列
            var = new double[new_data_X.Length];
            double[][] C = new double[I][];
            for (int i = 0; i < I; i++)
                C[i] = new double[I];
            double[,] C_inv;
            double[] c_vec = new double[I];
            double c_sca;

            for (int i = 0; i < I; i++)
                for (int j = 0; j < I; j++)
                    C[i][j] = GaussianKernel(data_X[i], data_X[j]) / alpha + 1 / beta;
            C_inv = jagToarray(C);
            C_inv = Mt.Inverse(jagToarray(C));  //Inner2zikeisiki

            for (int n = 0; n < new_data_X.Length; n++)
            {
                for (int i = 0; i < I; i++)
                    c_vec[i] = GaussianKernel(data_X[i], new_data_X[n]) / alpha;
                c_sca = GaussianKernel(new_data_X[n], new_data_X[n]) / alpha + 1 / beta;

                double[] y = new double[I];
                for (int i = 0; i < I; i++)
                    y[i] = data_y[i];


                exp[n] = Mt.Inner(C_inv, c_vec, y);
                var[n] = c_sca - Mt.Inner(C_inv, c_vec);

            }

        }
        static double GaussianKernel(double[] x, double[] x_dash)
        {
            double output;
            double norm2 = Mt.Norm2(Mt.Sub(x, x_dash));
            norm2 = norm2 * norm2;
            output = Math.Exp(-norm2 / (2 * 0.25));
            return output;
        }
        static double[,] jagToarray(double[][] input)
        {
            double[,] output = new double[input.Length, input[0].Length];
            for (int i = 0; i < input.Length; i++)
            {
                for (int j = 0; j < input[0].Length; j++)
                {
                    output[i, j] = input[i][j];
                }
            }
            return output;
        }
        static double[][] arrayTojag(double[,] input)
        {
            double[][] output = new double[input.GetLength(0)][];
            for (int i = 0; i < input.GetLength(0); i++)
            {
                output[i] = new double[input.GetLength(1)];
                for (int j = 0; j < input.GetLength(1); j++)
                {
                    output[i][j] = input[i, j];
                }
            }
            return output;
        }
    }
}
