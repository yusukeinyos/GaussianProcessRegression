using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace CsvFileIO
{
    public static class CsvFileIO
    {
        //----------------------------------------------------------------------------------------
        //行列ファイルの読み込み
        public static double[,] ReadData(string filename)
        {
            // ファイルを読み込むだけ
            // インスタンスを作成，パスをコンストラクタに渡す
            //StreamReader sr = new StreamReader("data\\MAX120.txt", Encoding.GetEncoding("Shift_JIS"));
            StreamReader sr = new StreamReader(filename, Encoding.GetEncoding("Shift_JIS"));
            string sLine = ""; // ←　一時格納用
            List<string> TextData = new List<string>();
            string[][] num;

            int I;
            int J;

            double[,] X;

            // ファイルから1行ずつ読み込む
            while (sLine != null)
            {
                // ファイルから1行ずつ読み込むReadToEndでもいい
                sLine = sr.ReadLine();
                if (sLine != null)
                    TextData.Add(sLine); // addメソッドで追記
            }
            sr.Close();


            // 行列数を数える
            int line_count = TextData.Count; // 行数を数えているだけ
            I = line_count;
            // 列数を数えるために適当な配列に格納して数えてみる
            string temp = (string)TextData[0]; // string型にキャスト
            // splitメソッドで文字列アレイにして数える
            string[] temp2 = temp.Split('\t');
            int col_count = temp2.Length; // 列数を数えているだけ
            J = col_count;
            // -----------------------------------

            // ArrayListから2次元配列に格納する．
            // 2次元配列の定義
            num = new string[line_count][];
            for (int i = 0; i < line_count; i++)
                num[i] = new string[col_count];

            int a = 0, b = 0;
            foreach (string sOutput in TextData)
            {
                // 一行ずつ読み込んで，各行をsplitメソッドで分ける
                string[] temp_line = sOutput.Split('\t');
                foreach (string value in temp_line)
                {
                    num[a][b] = Convert.ToString(value);
                    b++;
                }
                b = 0;
                a++;
            }

            X = new double[I, J];

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    X[i, j] = Double.Parse(num[i][j]);

            return X;
        }
        //----------------------------------------------------------------------------------------
        public static void WriteData(string filename, double[,] data)
        {
            int I = data.GetLength(0);
            int J = data.GetLength(1);
            string temp;

            StreamWriter swT = new StreamWriter(filename);
            for (int i = 0; i < I; i++)
            {
                temp = "";
                for (int j = 0; j < J; j++)
                {
                    //temp += T[i,k].ToString("f2");
                    temp += data[i, j].ToString();
                    if (j != J - 1) temp += (",");
                }
                swT.WriteLine(temp);
            }
            swT.Close();
        }
        //----------------------------------------------------------------------------------------
    }
}
