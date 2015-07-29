using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;  //プロジェクトの参照設定でSystem.Numericsを追加
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;
//プロジェクトのビルド設定でunsafe codeを許可

namespace System.Linq
{
    using InoueLab;

    #region Extension class
    public static partial class SystemLinqEx
    {
        #region Int32
        public static double GeometricalAverage<TSource>(this IEnumerable<TSource> source, Func<TSource, int> selector) { return source.Select(element => selector(element)).GeometricalAverage(); }
        public static int Product<TSource>(this IEnumerable<TSource> source, Func<TSource, int> selector) { return source.Select(element => selector(element)).Product(); }
        public static double GeometricalAverage(this IEnumerable<int> source) { var r = source._Product(); return Math.Pow((double)r.v0, 1.0 / r.v1); }
        public static int Product(this IEnumerable<int> source) { return source._Product().v0; }
        static V2<int, int> _Product(this IEnumerable<int> source)
        {
            int a = 1;
            int c = 0;
            foreach (var e in source) { checked { a *= e; } c++; }
            return New.V2(a, c);
        }
        #endregion

        #region Double
        public static double GeometricalAverage<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).GeometricalAverage(); }
        public static double Product<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).Product(); }
        public static double GeometricalAverage(this IEnumerable<double> source) { var r = source._Product(); return Math.Pow(r.v0, 1.0 / r.v1); }
        public static double Product(this IEnumerable<double> source) { return source._Product().v0; }
        static V2<double, int> _Product(this IEnumerable<double> source)
        {
            double a = 1;
            int c = 0;
            foreach (var e in source) { a *= e; c++; }
            return New.V2(a, c);
        }
        #endregion

        #region UInt32
        public static double Average<TSource>(this IEnumerable<TSource> source, Func<TSource, uint> selector) { return source.Select(element => selector(element)).Average(); }
        public static uint Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, uint> selector) { return source.Select(element => selector(element)).Sum(); }
        public static double Average(this IEnumerable<uint> source) { var r = source._Sum(); return (double)r.v0 / r.v1; }
        public static uint Sum(this IEnumerable<uint> source) { return source._Sum().v0; }
        static V2<uint, int> _Sum(this IEnumerable<uint> source)
        {
            var a = default(uint);
            int c = 0;
            foreach (var e in source) { checked { a += e; } c++; }
            return New.V2(a, c);
        }
        public static double GeometricalAverage<TSource>(this IEnumerable<TSource> source, Func<TSource, uint> selector) { return source.Select(element => selector(element)).GeometricalAverage(); }
        public static uint Product<TSource>(this IEnumerable<TSource> source, Func<TSource, uint> selector) { return source.Select(element => selector(element)).Product(); }
        public static double GeometricalAverage(this IEnumerable<uint> source) { var r = source._Product(); return Math.Pow((double)r.v0, 1.0 / r.v1); }
        public static uint Product(this IEnumerable<uint> source) { return source._Product().v0; }
        static V2<uint, int> _Product(this IEnumerable<uint> source)
        {
            uint a = 1;
            int c = 0;
            foreach (var e in source) { checked { a *= e; } c++; }
            return New.V2(a, c);
        }
        #endregion

        #region UInt64
        public static double Average<TSource>(this IEnumerable<TSource> source, Func<TSource, ulong> selector) { return source.Select(element => selector(element)).Average(); }
        public static ulong Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, ulong> selector) { return source.Select(element => selector(element)).Sum(); }
        public static double Average(this IEnumerable<ulong> source) { var r = source._Sum(); return (double)r.v0 / r.v1; }
        public static ulong Sum(this IEnumerable<ulong> source) { return source._Sum().v0; }
        static V2<ulong, int> _Sum(this IEnumerable<ulong> source)
        {
            var a = default(ulong);
            int c = 0;
            foreach (var e in source) { checked { a += e; } c++; }
            return New.V2(a, c);
        }
        public static double GeometricalAverage<TSource>(this IEnumerable<TSource> source, Func<TSource, ulong> selector) { return source.Select(element => selector(element)).GeometricalAverage(); }
        public static ulong Product<TSource>(this IEnumerable<TSource> source, Func<TSource, ulong> selector) { return source.Select(element => selector(element)).Product(); }
        public static double GeometricalAverage(this IEnumerable<ulong> source) { var r = source._Product(); return Math.Pow((double)r.v0, 1.0 / r.v1); }
        public static ulong Product(this IEnumerable<ulong> source) { return source._Product().v0; }
        static V2<ulong, int> _Product(this IEnumerable<ulong> source)
        {
            ulong a = 1;
            int c = 0;
            foreach (var e in source) { checked { a *= e; } c++; }
            return New.V2(a, c);
        }
        #endregion

        #region BigInteger
        public static double Average<TSource>(this IEnumerable<TSource> source, Func<TSource, BigInteger> selector) { return source.Select(element => selector(element)).Average(); }
        public static BigInteger Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, BigInteger> selector) { return source.Select(element => selector(element)).Sum(); }
        public static double Average(this IEnumerable<BigInteger> source) { var r = source._Sum(); return (double)r.v0 / r.v1; }
        public static BigInteger Sum(this IEnumerable<BigInteger> source) { return source._Sum().v0; }
        static V2<BigInteger, int> _Sum(this IEnumerable<BigInteger> source)
        {
            var a = default(BigInteger);
            int c = 0;
            foreach (var e in source) { a += e; c++; }
            return New.V2(a, c);
        }
        public static double GeometricalAverage<TSource>(this IEnumerable<TSource> source, Func<TSource, BigInteger> selector) { return source.Select(element => selector(element)).GeometricalAverage(); }
        public static BigInteger Product<TSource>(this IEnumerable<TSource> source, Func<TSource, BigInteger> selector) { return source.Select(element => selector(element)).Product(); }
        public static double GeometricalAverage(this IEnumerable<BigInteger> source) { var r = source._Product(); return Math.Pow((double)r.v0, 1.0 / r.v1); }
        public static BigInteger Product(this IEnumerable<BigInteger> source) { return source._Product().v0; }
        static V2<BigInteger, int> _Product(this IEnumerable<BigInteger> source)
        {
            BigInteger a = 1;
            int c = 0;
            foreach (var e in source) { a *= e; c++; }
            return New.V2(a, c);
        }
        #endregion

        #region TimeSpan
        public static TimeSpan Average<TSource>(this IEnumerable<TSource> source, Func<TSource, TimeSpan> selector) { return source.Select(element => selector(element)).Average(); }
        public static TimeSpan Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, TimeSpan> selector) { return source.Select(element => selector(element)).Sum(); }
        public static TimeSpan Average(this IEnumerable<TimeSpan> source) { var r = source._Sum(); return TimeSpan.FromMilliseconds(r.v0.TotalMilliseconds / r.v1); }
        public static TimeSpan Sum(this IEnumerable<TimeSpan> source) { return source._Sum().v0; }
        static V2<TimeSpan, int> _Sum(this IEnumerable<TimeSpan> source)
        {
            var a = default(TimeSpan);
            int c = 0;
            foreach (var e in source) { checked { a += e; } c++; }
            return New.V2(a, c);
        }
        #endregion

        #region Int2, Int3, Double2, Double3
        public static Double2 Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Int2> selector) { return source.Select(element => selector(element)).Average(); }
        public static Int2 Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, Int2> selector) { return source.Select(element => selector(element)).Sum(); }
        public static Double2 Average(this IEnumerable<Int2> source) { var r = source._Sum(); return (Double2)r.v0 / r.v1; }
        public static Int2 Sum(this IEnumerable<Int2> source) { return source._Sum().v0; }
        static V2<Int2, int> _Sum(this IEnumerable<Int2> source)
        {
            var a = default(Int2);
            int c = 0;
            foreach (var e in source) { a += e; c++; }
            return New.V2(a, c);
        }

        public static Double3 Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Int3> selector) { return source.Select(element => selector(element)).Average(); }
        public static Int3 Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, Int3> selector) { return source.Select(element => selector(element)).Sum(); }
        public static Double3 Average(this IEnumerable<Int3> source) { var r = source._Sum(); return (Double3)r.v0 / r.v1; }
        public static Int3 Sum(this IEnumerable<Int3> source) { return source._Sum().v0; }
        static V2<Int3, int> _Sum(this IEnumerable<Int3> source)
        {
            var a = default(Int3);
            int c = 0;
            foreach (var e in source) { a += e; c++; }
            return New.V2(a, c);
        }

        public static Double2 Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Double2> selector) { return source.Select(element => selector(element)).Average(); }
        public static Double2 Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, Double2> selector) { return source.Select(element => selector(element)).Sum(); }
        public static Double2 Average(this IEnumerable<Double2> source) { var r = source._Sum(); return r.v0 / r.v1; }
        public static Double2 Sum(this IEnumerable<Double2> source) { return source._Sum().v0; }
        static V2<Double2, int> _Sum(this IEnumerable<Double2> source)
        {
            var a = default(Double2);
            int c = 0;
            foreach (var e in source) { a += e; c++; }
            return New.V2(a, c);
        }

        public static Double3 Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Double3> selector) { return source.Select(element => selector(element)).Average(); }
        public static Double3 Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, Double3> selector) { return source.Select(element => selector(element)).Sum(); }
        public static Double3 Average(this IEnumerable<Double3> source) { var r = source._Sum(); return r.v0 / r.v1; }
        public static Double3 Sum(this IEnumerable<Double3> source) { return source._Sum().v0; }
        static V2<Double3, int> _Sum(this IEnumerable<Double3> source)
        {
            var a = default(Double3);
            int c = 0;
            foreach (var e in source) { a += e; c++; }
            return New.V2(a, c);
        }
        #endregion

        #region Double[]
        public static double[] Average<TSource>(this IEnumerable<TSource> source, Func<TSource, double[]> selector) { return source.Select(element => selector(element)).Average(); }
        public static double[] Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, double[]> selector) { return source.Select(element => selector(element)).Sum(); }
        public static double[] Average(this IEnumerable<double[]> source) { var r = source._Sum(); return r.v0.LetDiv(r.v1); }
        public static double[] Sum(this IEnumerable<double[]> source) { return source._Sum().v0; }
        static V2<double[], int> _Sum(this IEnumerable<double[]> source)
        {
            var a = default(double[]);
            int c = 0;
            foreach (var e in source) { if (c == 0) a = e.CloneX(); else a.LetAdd(e); c++; }
            return New.V2(a, c);
        }
        #endregion

        #region Double[,]
        public static double[,] Average<TSource>(this IEnumerable<TSource> source, Func<TSource, double[,]> selector) { return source.Select(element => selector(element)).Average(); }
        public static double[,] Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, double[,]> selector) { return source.Select(element => selector(element)).Sum(); }
        public static double[,] Average(this IEnumerable<double[,]> source) { var r = source._Sum(); return r.v0.LetDiv(r.v1); }
        public static double[,] Sum(this IEnumerable<double[,]> source) { return source._Sum().v0; }
        static V2<double[,], int> _Sum(this IEnumerable<double[,]> source)
        {
            var a = default(double[,]);
            int c = 0;
            foreach (var e in source) { if (c == 0) a = e.CloneX(); else a.LetAdd(e); c++; }
            return New.V2(a, c);
        }
        #endregion

        #region complex
        public static complex Average<TSource>(this IEnumerable<TSource> source, Func<TSource, complex> selector) { return source.Select(element => selector(element)).Average(); }
        public static complex Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, complex> selector) { return source.Select(element => selector(element)).Sum(); }
        public static complex Average(this IEnumerable<complex> source) { var r = source._Sum(); return r.v0 / r.v1; }
        public static complex Sum(this IEnumerable<complex> source) { return source._Sum().v0; }
        static V2<complex, int> _Sum(this IEnumerable<complex> source)
        {
            var a = default(complex);
            int c = 0;
            foreach (var e in source) { a += e; c++; }
            return New.V2(a, c);
        }
        public static complex GeometricalAverage<TSource>(this IEnumerable<TSource> source, Func<TSource, complex> selector) { return source.Select(element => selector(element)).GeometricalAverage(); }
        public static complex Product<TSource>(this IEnumerable<TSource> source, Func<TSource, complex> selector) { return source.Select(element => selector(element)).Product(); }
        public static complex GeometricalAverage(this IEnumerable<complex> source) { var r = source._Product(); return complex.Pow(r.v0, 1.0 / r.v1); }
        public static complex Product(this IEnumerable<complex> source) { return source._Product().v0; }
        static V2<complex, int> _Product(this IEnumerable<complex> source)
        {
            complex a = complex.One;
            int c = 0;
            foreach (var e in source) { a *= e; c++; }
            return New.V2(a, c);
        }
        #endregion

        #region complex[]
        public static complex[] Average<TSource>(this IEnumerable<TSource> source, Func<TSource, complex[]> selector) { return source.Select(element => selector(element)).Average(); }
        public static complex[] Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, complex[]> selector) { return source.Select(element => selector(element)).Sum(); }
        public static complex[] Average(this IEnumerable<complex[]> source) { var r = source._Sum(); return r.v0.LetDiv(r.v1); }
        public static complex[] Sum(this IEnumerable<complex[]> source) { return source._Sum().v0; }
        static V2<complex[], int> _Sum(this IEnumerable<complex[]> source)
        {
            var a = default(complex[]);
            int c = 0;
            foreach (var e in source) { if (c == 0) a = e.CloneX(); else a.LetAdd(e); c++; }
            return New.V2(a, c);
        }
        #endregion

        #region complex[,]
        public static complex[,] Average<TSource>(this IEnumerable<TSource> source, Func<TSource, complex[,]> selector) { return source.Select(element => selector(element)).Average(); }
        public static complex[,] Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, complex[,]> selector) { return source.Select(element => selector(element)).Sum(); }
        public static complex[,] Average(this IEnumerable<complex[,]> source) { var r = source._Sum(); return r.v0.LetDiv(r.v1); }
        public static complex[,] Sum(this IEnumerable<complex[,]> source) { return source._Sum().v0; }
        static V2<complex[,], int> _Sum(this IEnumerable<complex[,]> source)
        {
            var a = default(complex[,]);
            int c = 0;
            foreach (var e in source) { if (c == 0) a = e.CloneX(); else a.LetAdd(e); c++; }
            return New.V2(a, c);
        }
        #endregion
    }
    #endregion
}

namespace InoueLab
{
    #region Error class
    public static class ThrowException
    {
        public static void WriteLine(string message)
        {
            Console.Error.WriteLine(message);
        }
        public static void NotImplemented()
        {
            throw new NotImplementedException();
        }
        public static T NotImplemented<T>()
        {
            throw new NotImplementedException();
        }
        public static void Argument(string message)
        {
            throw new ArgumentException(message);
        }
        public static T Argument<T>(string message)
        {
            throw new ArgumentException(message);
        }
        public static void ArgumentNull(string paramName)
        {
            throw new ArgumentNullException(paramName);
        }
        public static void ArgumentOutOfRange(string paramName)
        {
            throw new ArgumentOutOfRangeException(paramName);
        }
        public static T ArgumentOutOfRange<T>(string paramName)
        {
            throw new ArgumentOutOfRangeException(paramName);
        }
        public static void InvalidOperation(string message)
        {
            throw new InvalidOperationException(message);
        }
        public static void NoElements()
        {
            throw new InvalidOperationException("no elements");
        }
        public static void SizeMismatch()
        {
            throw new InvalidOperationException("size mismatch");
        }
    }
    #endregion

    #region Concatenated structs
    [Serializable]
    public struct V2<T0, T1>
    {
        public T0 v0;
        public T1 v1;
        public V2(T0 v0, T1 v1)
        {
            this.v0 = v0;
            this.v1 = v1;
        }
        public override int GetHashCode() { return Ex.CombineHashCodes(v0.GetHashCode(), v1.GetHashCode()); }
        public override string ToString() { return string.Format("{0}, {1}", v0, v1); }
    }

    [Serializable]
    public struct V3<T0, T1, T2>
    {
        public T0 v0;
        public T1 v1;
        public T2 v2;
        public V3(T0 v0, T1 v1, T2 v2)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
        }
        public override int GetHashCode() { return Ex.CombineHashCodes(v0.GetHashCode(), v1.GetHashCode(), v2.GetHashCode()); }
        public override string ToString() { return string.Format("{0}, {1}, {2}", v0, v1, v2); }
    }

    [Serializable]
    public struct V4<T0, T1, T2, T3>
    {
        public T0 v0;
        public T1 v1;
        public T2 v2;
        public T3 v3;
        public V4(T0 v0, T1 v1, T2 v2, T3 v3)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
        }
        public override int GetHashCode() { return Ex.CombineHashCodes(v0.GetHashCode(), v1.GetHashCode(), v2.GetHashCode(), v3.GetHashCode()); }
        public override string ToString() { return string.Format("{0}, {1}, {2}, {3}", v0, v1, v2, v3); }
    }
    #endregion

    #region Fixed-size arrays
    [Serializable]
    public struct Array1<T> : IList<T>
    {
        public T v0;
        public Array1(T v0)
        {
            this.v0 = v0;
        }
        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            return -1;
        }
        public void Insert(int index, T item) { ThrowException.NotImplemented(); }
        public void RemoveAt(int index) { ThrowException.NotImplemented(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                }
                return ThrowException.Argument<T>("index");
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                }
                ThrowException.Argument("index");
            }
        }
        public void Add(T item) { ThrowException.NotImplemented(); }
        public void Clear() { ThrowException.NotImplemented(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex] = v0;
        }
        public int Count { get { return 1; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { return ThrowException.NotImplemented<bool>(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }

    [Serializable]
    public struct Array2<T> : IList<T>
    {
        public T v0, v1;
        public Array2(T v0, T v1)
        {
            this.v0 = v0;
            this.v1 = v1;
        }
        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            if (item.Equals(v1)) return 1;
            return -1;
        }
        public void Insert(int index, T item) { ThrowException.NotImplemented(); }
        public void RemoveAt(int index) { ThrowException.NotImplemented(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                }
                return ThrowException.Argument<T>("index");
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                    case 1: v1 = value; return;
                }
                ThrowException.Argument("index");
            }
        }
        public void Add(T item) { ThrowException.NotImplemented(); }
        public void Clear() { ThrowException.NotImplemented(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex + 0] = v0;
            array[arrayIndex + 1] = v1;
        }
        public int Count { get { return 2; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { return ThrowException.NotImplemented<bool>(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
            yield return v1;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }

    [Serializable]
    public struct Array3<T> : IList<T>
    {
        public T v0, v1, v2;
        public Array3(T v0, T v1, T v2)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
        }
        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            if (item.Equals(v1)) return 1;
            if (item.Equals(v2)) return 2;
            return -1;
        }
        public void Insert(int index, T item) { ThrowException.NotImplemented(); }
        public void RemoveAt(int index) { ThrowException.NotImplemented(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                    case 2: return v2;
                }
                return ThrowException.Argument<T>("index");
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                    case 1: v1 = value; return;
                    case 2: v2 = value; return;
                }
                ThrowException.Argument("index");
            }
        }
        public void Add(T item) { ThrowException.NotImplemented(); }
        public void Clear() { ThrowException.NotImplemented(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex + 0] = v0;
            array[arrayIndex + 1] = v1;
            array[arrayIndex + 2] = v2;
        }
        public int Count { get { return 2; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { return ThrowException.NotImplemented<bool>(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
            yield return v1;
            yield return v2;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }

    [Serializable]
    public struct Array4<T> : IList<T>
    {
        public T v0, v1, v2, v3;
        public Array4(T v0, T v1, T v2, T v3)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
        }
        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            if (item.Equals(v1)) return 1;
            if (item.Equals(v2)) return 2;
            if (item.Equals(v3)) return 3;
            return -1;
        }
        public void Insert(int index, T item) { ThrowException.NotImplemented(); }
        public void RemoveAt(int index) { ThrowException.NotImplemented(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                    case 2: return v2;
                    case 3: return v3;
                }
                return ThrowException.Argument<T>("index");
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                    case 1: v1 = value; return;
                    case 2: v2 = value; return;
                    case 3: v3 = value; return;
                }
                ThrowException.Argument("index");
            }
        }
        public void Add(T item) { ThrowException.NotImplemented(); }
        public void Clear() { ThrowException.NotImplemented(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex + 0] = v0;
            array[arrayIndex + 1] = v1;
            array[arrayIndex + 2] = v2;
            array[arrayIndex + 3] = v3;
        }
        public int Count { get { return 2; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { return ThrowException.NotImplemented<bool>(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
            yield return v1;
            yield return v2;
            yield return v3;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }
    #endregion

    #region Int2, Int3, Double2, Double3
    [Serializable]
    public partial struct Int2 : IEquatable<Int2>, IComparable<Int2>, IFormattable
    {
        public int v0, v1;
        public int X { get { return v0; } set { v0 = value; } }
        public int Y { get { return v1; } set { v1 = value; } }
        public int this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                }
                return ThrowException.ArgumentOutOfRange<int>("index");
            }
        }
        public Int2(int v0, int v1) { this.v0 = v0; this.v1 = v1; }

        public static explicit operator Int3(Int2 value) { return new Int3(value.v0, value.v1, 0); }
        public static explicit operator Double2(Int2 value) { return new Double2(value.v0, value.v1); }
        public static explicit operator Double3(Int2 value) { return new Double3(value.v0, value.v1, 0); }
        public static bool operator ==(Int2 left, Int2 right) { return ((left.v0 - right.v0) | (left.v1 - right.v1)) == 0; }
        public static bool operator !=(Int2 left, Int2 right) { return ((left.v0 - right.v0) | (left.v1 - right.v1)) != 0; }
        public static bool operator <=(Int2 left, Int2 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Int2 left, Int2 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Int2 left, Int2 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Int2 left, Int2 right) { return left.CompareTo(right) > 0; }
        public static Int2 operator ~(Int2 value) { return new Int2(~value.v0, ~value.v1); }
        public static Int2 operator +(Int2 value) { return new Int2(+value.v0, +value.v1); }
        public static Int2 operator -(Int2 value) { return new Int2(-value.v0, -value.v1); }
        public static Int2 operator +(Int2 left, Int2 right) { return new Int2(left.v0 + right.v0, left.v1 + right.v1); }
        public static Int2 operator -(Int2 left, Int2 right) { return new Int2(left.v0 - right.v0, left.v1 - right.v1); }
        public static Int2 operator *(Int2 left, Int2 right) { return new Int2(left.v0 * right.v0, left.v1 * right.v1); }
        public static Int2 operator /(Int2 left, Int2 right) { return new Int2(left.v0 / right.v0, left.v1 / right.v1); }
        public static Int2 operator %(Int2 left, Int2 right) { return new Int2(left.v0 % right.v0, left.v1 % right.v1); }
        public static Int2 operator |(Int2 left, Int2 right) { return new Int2(left.v0 | right.v0, left.v1 | right.v1); }
        public static Int2 operator &(Int2 left, Int2 right) { return new Int2(left.v0 & right.v0, left.v1 & right.v1); }
        public static Int2 operator ^(Int2 left, Int2 right) { return new Int2(left.v0 ^ right.v0, left.v1 ^ right.v1); }
        public static Int2 operator |(Int2 left, int right) { return new Int2(left.v0 | right, left.v1 | right); }
        public static Int2 operator &(Int2 left, int right) { return new Int2(left.v0 & right, left.v1 & right); }
        public static Int2 operator ^(Int2 left, int right) { return new Int2(left.v0 ^ right, left.v1 ^ right); }
        public static Int2 operator +(Int2 left, int right) { return new Int2(left.v0 + right, left.v1 + right); }
        public static Int2 operator -(Int2 left, int right) { return new Int2(left.v0 - right, left.v1 - right); }
        public static Int2 operator *(Int2 left, int right) { return new Int2(left.v0 * right, left.v1 * right); }
        public static Int2 operator /(Int2 left, int right) { return new Int2(left.v0 / right, left.v1 / right); }
        public static Int2 operator %(Int2 left, int right) { return new Int2(left.v0 % right, left.v1 % right); }
        public override bool Equals(object obj) { return (obj is Int2) && (this == (Int2)obj); }
        public override int GetHashCode() { return Ex.CombineHashCodes(v0.GetHashCode(), v1.GetHashCode()); }
        public override string ToString() { return string.Format("({0}, {1})", v0, v1); }
        public string ToString(string format) { return string.Format("({0}, {1})", v0.ToString(format), v1.ToString(format)); }
        public string ToString(IFormatProvider provider) { return string.Format(provider, "({0}, {1})", v0, v1); }
        public string ToString(string format, IFormatProvider provider) { return string.Format(provider, "({0}, {1})", v0.ToString(format, provider), v1.ToString(format, provider)); }
        public bool Equals(Int2 other) { return ((v0 - other.v0) | (v1 - other.v1)) == 0; }
        public int CompareTo(Int2 other) { return v1 != other.v1 ? v1 - other.v1 : v0 - other.v0; }
    }

    [Serializable]
    public partial struct Int3 : IEquatable<Int3>, IComparable<Int3>, IFormattable
    {
        public int v0, v1, v2;
        public int X { get { return v0; } set { v0 = value; } }
        public int Y { get { return v1; } set { v1 = value; } }
        public int Z { get { return v2; } set { v2 = value; } }
        public int this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                    case 2: return v2;
                }
                return ThrowException.ArgumentOutOfRange<int>("index");
            }
        }
        public Int3(int v0, int v1, int v2) { this.v0 = v0; this.v1 = v1; this.v2 = v2; }

        public static explicit operator Int2(Int3 value) { return new Int2(value.v0, value.v1); }
        public static explicit operator Double2(Int3 value) { return new Double2(value.v0, value.v1); }
        public static explicit operator Double3(Int3 value) { return new Double3(value.v0, value.v1, value.v2); }
        public static bool operator ==(Int3 left, Int3 right) { return ((left.v0 - right.v0) | (left.v1 - right.v1) | (left.v2 - right.v2)) == 0; }
        public static bool operator !=(Int3 left, Int3 right) { return ((left.v0 - right.v0) | (left.v1 - right.v1) | (left.v2 - right.v2)) != 0; }
        public static bool operator <=(Int3 left, Int3 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Int3 left, Int3 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Int3 left, Int3 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Int3 left, Int3 right) { return left.CompareTo(right) > 0; }
        public static Int3 operator ~(Int3 value) { return new Int3(~value.v0, ~value.v1, ~value.v2); }
        public static Int3 operator +(Int3 value) { return new Int3(+value.v0, +value.v1, +value.v2); }
        public static Int3 operator -(Int3 value) { return new Int3(-value.v0, -value.v1, -value.v2); }
        public static Int3 operator +(Int3 left, Int3 right) { return new Int3(left.v0 + right.v0, left.v1 + right.v1, left.v2 + right.v2); }
        public static Int3 operator -(Int3 left, Int3 right) { return new Int3(left.v0 - right.v0, left.v1 - right.v1, left.v2 - right.v2); }
        public static Int3 operator *(Int3 left, Int3 right) { return new Int3(left.v0 * right.v0, left.v1 * right.v1, left.v2 * right.v2); }
        public static Int3 operator /(Int3 left, Int3 right) { return new Int3(left.v0 / right.v0, left.v1 / right.v1, left.v2 / right.v2); }
        public static Int3 operator %(Int3 left, Int3 right) { return new Int3(left.v0 % right.v0, left.v1 % right.v1, left.v2 % right.v2); }
        public static Int3 operator |(Int3 left, Int3 right) { return new Int3(left.v0 | right.v0, left.v1 | right.v1, left.v2 | right.v2); }
        public static Int3 operator &(Int3 left, Int3 right) { return new Int3(left.v0 & right.v0, left.v1 & right.v1, left.v2 & right.v2); }
        public static Int3 operator ^(Int3 left, Int3 right) { return new Int3(left.v0 ^ right.v0, left.v1 ^ right.v1, left.v2 ^ right.v2); }
        public static Int3 operator |(Int3 left, int right) { return new Int3(left.v0 | right, left.v1 | right, left.v2 | right); }
        public static Int3 operator &(Int3 left, int right) { return new Int3(left.v0 & right, left.v1 & right, left.v2 & right); }
        public static Int3 operator ^(Int3 left, int right) { return new Int3(left.v0 ^ right, left.v1 ^ right, left.v2 ^ right); }
        public static Int3 operator +(Int3 left, int right) { return new Int3(left.v0 + right, left.v1 + right, left.v2 + right); }
        public static Int3 operator -(Int3 left, int right) { return new Int3(left.v0 - right, left.v1 - right, left.v2 - right); }
        public static Int3 operator *(Int3 left, int right) { return new Int3(left.v0 * right, left.v1 * right, left.v2 * right); }
        public static Int3 operator /(Int3 left, int right) { return new Int3(left.v0 / right, left.v1 / right, left.v2 / right); }
        public static Int3 operator %(Int3 left, int right) { return new Int3(left.v0 % right, left.v1 % right, left.v2 % right); }
        public override bool Equals(object obj) { return (obj is Int3) && (this == (Int3)obj); }
        public override int GetHashCode() { return Ex.CombineHashCodes(v0.GetHashCode(), v1.GetHashCode(), v2.GetHashCode()); }
        public override string ToString() { return string.Format("({0}, {1}, {2})", v0, v1, v2); }
        public string ToString(string format) { return string.Format("({0}, {1}, {2})", v0.ToString(format), v1.ToString(format), v2.ToString(format)); }
        public string ToString(IFormatProvider provider) { return string.Format(provider, "({0}, {1}, {2})", v0, v1, v2); }
        public string ToString(string format, IFormatProvider provider) { return string.Format(provider, "({0}, {1}, {2})", v0.ToString(format, provider), v1.ToString(format, provider), v2.ToString(format, provider)); }
        public bool Equals(Int3 other) { return ((v0 - other.v0) | (v1 - other.v1) | (v2 - other.v2)) == 0; }
        public int CompareTo(Int3 other) { return v2 != other.v2 ? v2 - other.v2 : (v1 != other.v1 ? v1 - other.v1 : v0 - other.v0); }
    }

    [Serializable]
    public partial struct Double2 : IEquatable<Double2>, IComparable<Double2>, IFormattable
    {
        public static readonly Double2 Zero = new Double2(0, 0);
        public static readonly Double2 One = new Double2(1, 1);
        public double v0, v1;
        public double X { get { return v0; } set { v0 = value; } }
        public double Y { get { return v1; } set { v1 = value; } }
        public double this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                }
                return ThrowException.ArgumentOutOfRange<double>("index");
            }
        }
        public Double2(double v0, double v1) { this.v0 = v0; this.v1 = v1; }

        public static explicit operator Int2(Double2 value) { return new Int2((int)value.v0, (int)value.v1); }
        public static explicit operator Int3(Double2 value) { return new Int3((int)value.v0, (int)value.v1, 0); }
        public static explicit operator Double3(Double2 value) { return new Double3(value.v0, value.v1, 0); }
        public static bool operator ==(Double2 left, Double2 right) { return (left.v0 == right.v0) && (left.v1 == right.v1); }
        public static bool operator !=(Double2 left, Double2 right) { return (left.v0 != right.v0) || (left.v1 != right.v1); }
        public static bool operator <=(Double2 left, Double2 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Double2 left, Double2 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Double2 left, Double2 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Double2 left, Double2 right) { return left.CompareTo(right) > 0; }
        public static Double2 operator +(Double2 value) { return new Double2(+value.v0, +value.v1); }
        public static Double2 operator -(Double2 value) { return new Double2(-value.v0, -value.v1); }
        public static Double2 operator +(Double2 left, Double2 right) { return new Double2(left.v0 + right.v0, left.v1 + right.v1); }
        public static Double2 operator -(Double2 left, Double2 right) { return new Double2(left.v0 - right.v0, left.v1 - right.v1); }
        public static Double2 operator *(Double2 left, Double2 right) { return new Double2(left.v0 * right.v0, left.v1 * right.v1); }
        public static Double2 operator /(Double2 left, Double2 right) { return new Double2(left.v0 / right.v0, left.v1 / right.v1); }
        public static Double2 operator +(Double2 left, double right) { return new Double2(left.v0 + right, left.v1 + right); }
        public static Double2 operator -(Double2 left, double right) { return new Double2(left.v0 - right, left.v1 - right); }
        public static Double2 operator *(Double2 left, double right) { return new Double2(left.v0 * right, left.v1 * right); }
        public static Double2 operator /(Double2 left, double right) { return new Double2(left.v0 / right, left.v1 / right); }
        public override bool Equals(object obj) { return (obj is Double2) && (this == (Double2)obj); }
        public override int GetHashCode() { return Ex.CombineHashCodes(v0.GetHashCode(), v1.GetHashCode()); }
        public override string ToString() { return string.Format("({0}, {1})", v0, v1); }
        public string ToString(string format) { return string.Format("({0}, {1})", v0.ToString(format), v1.ToString(format)); }
        public string ToString(IFormatProvider provider) { return string.Format(provider, "({0}, {1})", v0, v1); }
        public string ToString(string format, IFormatProvider provider) { return string.Format(provider, "({0}, {1})", v0.ToString(format, provider), v1.ToString(format, provider)); }
        public bool Equals(Double2 other) { return (v0 == other.v0) && (v1 == other.v1); }
        public int CompareTo(Double2 other) { return v1 != other.v1 ? (v1 < other.v1 ? -1 : 1) : (v0 == other.v0 ? 0 : v0 < other.v0 ? -1 : 1); }
        public static bool IsNaN(Double2 value) { return double.IsNaN(value.v0) || double.IsNaN(value.v1); }
    }

    [Serializable]
    public partial struct Double3 : IEquatable<Double3>, IComparable<Double3>, IFormattable
    {
        public static readonly Double3 Zero = new Double3(0, 0, 0);
        public static readonly Double3 One = new Double3(1, 1, 1);
        public double v0, v1, v2;
        public double X { get { return v0; } set { v0 = value; } }
        public double Y { get { return v1; } set { v1 = value; } }
        public double Z { get { return v2; } set { v2 = value; } }
        public double this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                    case 2: return v2;
                }
                return ThrowException.ArgumentOutOfRange<double>("index");
            }
        }
        public Double3(double v0, double v1, double v2) { this.v0 = v0; this.v1 = v1; this.v2 = v2; }

        public static explicit operator Int2(Double3 value) { return new Int2((int)value.v0, (int)value.v1); }
        public static explicit operator Int3(Double3 value) { return new Int3((int)value.v0, (int)value.v1, (int)value.v2); }
        public static explicit operator Double2(Double3 value) { return new Double2(value.v0, value.v1); }
        public static bool operator ==(Double3 left, Double3 right) { return (left.v0 == right.v0) && (left.v1 == right.v1) && (left.v2 == right.v2); }
        public static bool operator !=(Double3 left, Double3 right) { return (left.v0 != right.v0) || (left.v1 != right.v1) || (left.v2 != right.v2); }
        public static bool operator <=(Double3 left, Double3 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Double3 left, Double3 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Double3 left, Double3 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Double3 left, Double3 right) { return left.CompareTo(right) > 0; }
        public static Double3 operator +(Double3 value) { return new Double3(+value.v0, +value.v1, +value.v2); }
        public static Double3 operator -(Double3 value) { return new Double3(-value.v0, -value.v1, -value.v2); }
        public static Double3 operator +(Double3 left, Double3 right) { return new Double3(left.v0 + right.v0, left.v1 + right.v1, left.v2 + right.v2); }
        public static Double3 operator -(Double3 left, Double3 right) { return new Double3(left.v0 - right.v0, left.v1 - right.v1, left.v2 - right.v2); }
        public static Double3 operator *(Double3 left, Double3 right) { return new Double3(left.v0 * right.v0, left.v1 * right.v1, left.v2 * right.v2); }
        public static Double3 operator /(Double3 left, Double3 right) { return new Double3(left.v0 / right.v0, left.v1 / right.v1, left.v2 / right.v2); }
        public static Double3 operator +(Double3 left, double right) { return new Double3(left.v0 + right, left.v1 + right, left.v2 + right); }
        public static Double3 operator -(Double3 left, double right) { return new Double3(left.v0 - right, left.v1 - right, left.v2 - right); }
        public static Double3 operator *(Double3 left, double right) { return new Double3(left.v0 * right, left.v1 * right, left.v2 * right); }
        public static Double3 operator /(Double3 left, double right) { return new Double3(left.v0 / right, left.v1 / right, left.v2 / right); }
        public override bool Equals(object obj) { return (obj is Double3) && (this == (Double3)obj); }
        public override int GetHashCode() { return Ex.CombineHashCodes(v0.GetHashCode(), v1.GetHashCode(), v2.GetHashCode()); }
        public override string ToString() { return string.Format("({0}, {1}, {2})", v0, v1, v2); }
        public string ToString(string format) { return string.Format("({0}, {1}, {2})", v0.ToString(format), v1.ToString(format), v2.ToString(format)); }
        public string ToString(IFormatProvider provider) { return string.Format(provider, "({0}, {1}, {2})", v0, v1, v2); }
        public string ToString(string format, IFormatProvider provider) { return string.Format(provider, "({0}, {1}, {2})", v0.ToString(format, provider), v1.ToString(format, provider), v2.ToString(format, provider)); }
        public bool Equals(Double3 other) { return (v0 == other.v0) && (v1 == other.v1) && (v2 == other.v2); }
        public int CompareTo(Double3 other) { return v2 != other.v2 ? (v2 < other.v1 ? -1 : 1) : (v1 != other.v1 ? (v1 < other.v1 ? -1 : 1) : (v0 == other.v0 ? 0 : v0 < other.v0 ? -1 : 1)); }
        public static bool IsNaN(Double3 value) { return double.IsNaN(value.v0) || double.IsNaN(value.v1) || double.IsNaN(value.v2); }
    }
    #endregion

    #region Ints class
    class Ints : IComparable, IEquatable<Ints>, IComparable<Ints>, IList<int>
    {
        readonly int[] _Items;
        public Ints(params int[] items)
        {
            this._Items = items.CloneX();
        }
        public Ints(IEnumerable<int> collections)
        {
            this._Items = collections.ToArray();
        }
        public int Length
        {
            get { return _Items.Length; }
        }
        public static bool operator ==(Ints left, Ints right) { return Equals(left, right); }
        public static bool operator !=(Ints left, Ints right) { return !Equals(left, right); }
        public static bool operator <=(Ints left, Ints right) { return Compare(left, right) <= 0; }
        public static bool operator >=(Ints left, Ints right) { return Compare(left, right) >= 0; }
        public static bool operator <(Ints left, Ints right) { return Compare(left, right) < 0; }
        public static bool operator >(Ints left, Ints right) { return Compare(left, right) > 0; }
        public static Ints operator +(Ints x, Ints y) { return new Ints(x._Items.Concat(y._Items)); }
        public static int Compare(Ints left, Ints right)
        {
            if ((object)left == null) return ((object)right == null) ? 0 : -1;
            return left.CompareTo(right);
        }
        public static bool Equals(Ints left, Ints right)
        {
            if ((object)left == null) return (object)right == null;
            return left.Equals(right);
        }
        public override bool Equals(object obj) { return (obj is Ints) && Equals((Ints)obj); }
        public override int GetHashCode()
        {
            int h = 0;
            for (int i = 0; i < _Items.Length; i++)
                h = h * 3 - 13 + _Items[i];
            return h;
        }
        public override string ToString() { return _Items.ToStringFormat("", ", "); }

        #region IComparable メンバー
        public int CompareTo(object obj)
        {
            var other = obj as Ints;
            if (other == null) ThrowException.Argument("obj");
            return CompareTo(other);
        }
        #endregion
        #region IEquatable<Ints> メンバ
        public bool Equals(Ints other)
        {
            if ((object)other == null) return false;
            if (_Items.Length != other._Items.Length) return false;
            for (int i = 0; i < _Items.Length; i++)
                if (_Items[i] != other._Items[i]) return false;
            return true;
        }
        #endregion
        #region IComparable<Ints> メンバ
        public int CompareTo(Ints other)
        {
            if ((object)other == null) return 1;
            int l = Math.Min(_Items.Length, other._Items.Length);
            for (int i = 0; i < l; i++)
                if (_Items[i] != other._Items[i]) return _Items[i] - other._Items[i];
            return _Items.Length - other._Items.Length;
        }
        #endregion
        #region IList<int> メンバー
        public int IndexOf(int item) { return _Items.IndexOf(item); }
        public void Insert(int index, int item) { ThrowException.NotImplemented(); }
        public void RemoveAt(int index) { ThrowException.NotImplemented(); }
        public int this[int index]
        {
            get { return _Items[index]; }
            set { ThrowException.NotImplemented(); }
        }
        #endregion
        #region ICollection<int> メンバー
        public void Add(int item) { ThrowException.NotImplemented(); }
        public void Clear() { ThrowException.NotImplemented(); }
        public bool Contains(int item) { return _Items.Contains(item); }
        public void CopyTo(int[] array, int arrayIndex) { _Items.CopyTo(array, arrayIndex); }
        public int Count
        {
            get { return _Items.Length; }
        }
        public bool IsReadOnly
        {
            get { return true; }
        }
        public bool Remove(int item) { return ThrowException.NotImplemented<bool>(); }
        #endregion
        #region IEnumerable<int> メンバー
        public IEnumerator<int> GetEnumerator()
        {
            for (int i = 0; i < _Items.Length; i++)
                yield return _Items[i];
        }
        #endregion
        #region IEnumerable メンバー
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
        #endregion
    }
    #endregion

    #region complex
    // System.Numerics.Complex
    [Serializable]
    public struct complex : IEquatable<complex>, IFormattable
    {
        public double Re;
        public double Im;
        public double Real { get { return Re; } }
        public double Imaginary { get { return Im; } }
        public double Magnitude { get { return Mt.Norm2(Re, Im); } }
        public double Phase { get { return Math.Atan2(Im, Re); } }
        public static readonly complex Zero = new complex(0, 0);
        public static readonly complex One = new complex(1, 0);
        public static readonly complex ImaginaryOne = new complex(0, 1);
        public complex(double real, double imaginary) { this.Re = real; this.Im = imaginary; }

        public static explicit operator complex(short value) { return new complex(value, 0); }
        public static explicit operator complex(int value) { return new complex(value, 0); }
        public static explicit operator complex(long value) { return new complex(value, 0); }
        public static explicit operator complex(ushort value) { return new complex(value, 0); }
        public static explicit operator complex(uint value) { return new complex(value, 0); }
        public static explicit operator complex(ulong value) { return new complex(value, 0); }
        public static explicit operator complex(sbyte value) { return new complex(value, 0); }
        public static explicit operator complex(byte value) { return new complex(value, 0); }
        public static explicit operator complex(float value) { return new complex(value, 0); }
        public static explicit operator complex(double value) { return new complex(value, 0); }
        public static explicit operator complex(decimal value) { return new complex((double)value, 0); }

        public static complex operator +(complex value) { return new complex(+value.Re, +value.Im); }
        public static complex operator -(complex value) { return new complex(-value.Re, -value.Im); }
        public static complex operator +(complex left, complex right) { return new complex(left.Re + right.Re, left.Im + right.Im); }
        public static complex operator -(complex left, complex right) { return new complex(left.Re - right.Re, left.Im - right.Im); }
        public static complex operator *(complex left, complex right) { return new complex(left.Re * right.Re - left.Im * right.Im, left.Re * right.Im + left.Im * right.Re); }
        public static complex operator /(complex left, complex right)
        {
            if (Math.Abs(right.Re) > Math.Abs(right.Im))
            {
                var a = right.Im / right.Re;
                var b = right.Re + right.Im * a;
                return new complex((left.Re + left.Im * a) / b, (left.Im - left.Re * a) / b);
            }
            else
            {
                var a = right.Re / right.Im;
                var b = right.Im + right.Re * a;
                return new complex((left.Im + left.Re * a) / b, (-left.Re + left.Im * a) / b);
            }
        }
        public static complex operator +(complex left, double right) { return new complex(left.Re + right, left.Im); }
        public static complex operator -(complex left, double right) { return new complex(left.Re - right, left.Im); }
        public static complex operator *(complex left, double right) { return new complex(left.Re * right, left.Im * right); }
        public static complex operator /(complex left, double right) { return new complex(left.Re / right, left.Im / right); }
        public static complex operator +(double left, complex right) { return new complex(left + right.Re, right.Im); }
        public static complex operator -(double left, complex right) { return new complex(left - right.Re, -right.Im); }
        public static complex operator *(double left, complex right) { return new complex(left * right.Re, left * right.Im); }
        public static complex operator /(double left, complex right)
        {
            if (Math.Abs(right.Re) > Math.Abs(right.Im))
            {
                var a = right.Im / right.Re;
                var b = right.Re + right.Im * a;
                return new complex(left / b, -left * a / b);
            }
            else
            {
                var a = right.Re / right.Im;
                var b = right.Im + right.Re * a;
                return new complex(left * a / b, -left / b);
            }
        }
        public static bool operator ==(complex left, complex right) { return left.Re == right.Re && left.Im == right.Im; }
        public static bool operator !=(complex left, complex right) { return left.Re != right.Re || left.Im != right.Im; }
        public static bool operator ==(complex left, double right) { return left.Re == right && left.Im == 0; }
        public static bool operator !=(complex left, double right) { return left.Re != right || left.Im != 0; }
        public static bool operator ==(double left, complex right) { return left == right.Re && 0 == right.Im; }
        public static bool operator !=(double left, complex right) { return left != right.Re || 0 != right.Im; }
        public override bool Equals(object obj) { return (obj is complex) ? this == (complex)obj : false; }
        public bool Equals(complex value) { return Re == value.Re && Im == value.Im; }
        public override string ToString() { return string.Format("({0}, {1})", Re, Im); }
        public string ToString(string format) { return string.Format("({0}, {1})", Re.ToString(format), Im.ToString(format)); }
        public string ToString(IFormatProvider provider) { return string.Format(provider, "({0}, {1})", Re, Im); }
        public string ToString(string format, IFormatProvider provider) { return string.Format(provider, "({0}, {1})", Re.ToString(format, provider), Im.ToString(format, provider)); }
        public override int GetHashCode()
        {
            return Re.GetHashCode() % 99999997 ^ Im.GetHashCode();
        }

        public static complex FromPolarCoordinates(double magnitude, double phase) { return new complex(magnitude * Math.Cos(phase), magnitude * Math.Sin(phase)); }
        public static double Abs(complex value) { return Mt.Norm2(value.Re, value.Im); }
        public static complex Conjugate(complex value) { return new complex(value.Re, -value.Im); }
        public static complex Reciprocal(complex value) { return (value.Re == 0 && value.Im == 0) ? value : 1 / value; }
        public static complex Cos(complex value) { return new complex(Math.Cos(value.Re) * Math.Cosh(value.Im), -Math.Sin(value.Re) * Math.Sinh(value.Im)); }
        public static complex Sin(complex value) { return new complex(Math.Sin(value.Re) * Math.Cosh(value.Im), Math.Cos(value.Re) * Math.Sinh(value.Im)); }
        public static complex Tan(complex value) { return Sin(value) / Cos(value); }
        public static complex Cosh(complex value) { return new complex(Math.Cosh(value.Re) * Math.Cos(value.Im), Math.Sinh(value.Re) * Math.Sin(value.Im)); }
        public static complex Sinh(complex value) { return new complex(Math.Sinh(value.Re) * Math.Cos(value.Im), Math.Cosh(value.Re) * Math.Sin(value.Im)); }
        public static complex Tanh(complex value) { return Sinh(value) / Cosh(value); }
        public static complex Acos(complex value)
        {
            //return -ImaginaryOne * Log(value + ImaginaryOne * Sqrt(One - value * value));
            var a = 1 - value.Re * value.Re + value.Im * value.Im;
            var b = -2 * value.Re * value.Im;
            var c = Math.Sqrt(Mt.Norm2(b, a));
            var d = Math.Atan2(b, a) / 2;
            var e = value.Re - Math.Sin(d) * c;
            var f = value.Im + Math.Cos(d) * c;
            return new complex(Math.Atan2(f, e), -Math.Log(Mt.Norm2(f, e)));
        }
        public static complex Asin(complex value)
        {
            //return -ImaginaryOne * Log(ImaginaryOne * value + Sqrt(One - value * value));
            var a = 1 - value.Re * value.Re + value.Im * value.Im;
            var b = -2 * value.Re * value.Im;
            var c = Math.Sqrt(Mt.Norm2(b, a));
            var d = Math.Atan2(b, a) / 2;
            var e = -value.Im + Math.Cos(d) * c;
            var f = +value.Re + Math.Sin(d) * c;
            return new complex(Math.Atan2(f, e), -Math.Log(Mt.Norm2(f, e)));
        }
        public static complex Atan(complex value)
        {
            //return new complex(0, 0.5) * (Log(One - ImaginaryOne * value) - Log(One + ImaginaryOne * value));
            return new complex(
                0.5 * (Math.Atan2(value.Re, 1 + value.Im) + Math.Atan2(value.Re, 1 - value.Im)),
                0.5 * Math.Log(Mt.Norm2(value.Re, 1 + value.Im) / Mt.Norm2(value.Re, 1 - value.Im))
            );
        }
        public static complex Log(complex value) { return new complex(Math.Log(Abs(value)), Math.Atan2(value.Im, value.Re)); }
        public static complex Log(complex value, double baseValue) { return Log(value) / Math.Log(baseValue); }
        public static complex Log10(complex value) { return Log(value) * 0.43429448190325; }
        public static complex Exp(complex value) { var f = Math.Exp(value.Re); return new complex(Math.Cos(value.Im) * f, Math.Sin(value.Im) * f); }
        public static complex Sqrt(complex value) { return FromPolarCoordinates(Math.Sqrt(value.Magnitude), value.Phase / 2); }
        public static complex Pow(complex value, complex power)
        {
            if (power.Im == 0) return Pow(value, power.Re);
            if (value == 0) return value;
            //return Exp(Log(value) * power);
            var a = Math.Log(Abs(value));
            var b = Math.Atan2(value.Im, value.Re);
            var c = Math.Exp(a * power.Re - b * power.Im);
            var d = /*******/a * power.Im + b * power.Re;
            return new complex(Math.Cos(d) * c, Math.Sin(d) * c);
        }
        public static complex Pow(complex value, double power)
        {
            if (power == 0) return One;
            if (power == 1 || value == 0) return value;
            return FromPolarCoordinates(Math.Pow(value.Magnitude, power), value.Phase * power);
        }
    }
    #endregion

    #region StopWatch
    public class StopWatch
    {
        DateTime RegisteredTime;
        TimeSpan Duration;
        bool Running;
        public TimeSpan Elapsed { get { return Running ? Duration + (DateTime.Now - RegisteredTime) : Duration; } }
        public bool IsRunning { get { return Running; } }
        public StopWatch() { Restart(); }
        public void Restart() { Reset(); Start(); }
        public void Reset() { Running = false; Duration = TimeSpan.Zero; }
        public void Start() { RegisteredTime = DateTime.Now; Running = true; }
        public TimeSpan Stop() { Running = false; return Duration += DateTime.Now - RegisteredTime; }
        public override string ToString() { return Elapsed.ToString(); }
    }
    #endregion

    #region RandomMT
    // MersenneTwister.dSFMT2.1
    [Serializable]
    public class RandomMT
    {
        const int N = 624, M = 397;
        const uint UPPER_MASK = 0x80000000u;
        const uint LOWER_MASK = 0x7fffffffu;
        const uint MATRIX_A = 0x9908b0dfu;

        uint BufferIndex;
        uint[] Buffer = new uint[N];

        public RandomMT() { Init((uint)DateTime.Now.ToBinary()); }
        public RandomMT(uint seed) { Init(seed); }
        public RandomMT(int seed) { Init((uint)seed); }

        public void Init(uint seed)
        {
            Buffer[0] = seed;
            for (int i = 1; i < N; i++)
                Buffer[i] = 1812433253u * (Buffer[i - 1] ^ (Buffer[i - 1] >> 30)) + (uint)i;
            BufferIndex = N;
        }
        public void Init(uint[] key)
        {
            if (key == null) ThrowException.ArgumentNull("key");
            Init(19650218u);
            uint i = 1, j = 0;
            for (int k = Math.Max(N, key.Length); k > 0; k--)
            {
                Buffer[i] = (Buffer[i] ^ ((Buffer[i - 1] ^ (Buffer[i - 1] >> 30)) * 1664525u)) + key[j] + j; // non linear
                if (++i >= N) { i = 1; Buffer[0] = Buffer[N - 1]; }
                if (++j >= key.Length) j = 0;
            }
            for (int k = N - 1; k > 0; k--)
            {
                Buffer[i] = (Buffer[i] ^ ((Buffer[i - 1] ^ (Buffer[i - 1] >> 30)) * 1566083941u)) - i; // non linear
                if (++i >= N) { i = 1; Buffer[0] = Buffer[N - 1]; }
            }
            Buffer[0] = 0x80000000u; // MSB is 1; assuring non-zero initial array
        }

        // 32-bit uint [0, 0xffffffff]
        public uint UInt32()
        {
            if (BufferIndex >= N) { BufferIndex = 0; UInt32_(); } // generate N words at one time            
            uint y = Buffer[BufferIndex++];
            // tempering
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680u;
            y ^= (y << 15) & 0xefc60000u;
            y ^= (y >> 18);
            return y;
        }
        void UInt32_()
        {
            for (int kk = 0, k1 = 1, kM = M; kk < N; kk++, k1++, kM++)
            {
                if (k1 == N) k1 = 0;
                if (kM == N) kM = 0;
                uint y = (Buffer[kk] & UPPER_MASK) | (Buffer[k1] & LOWER_MASK);
                Buffer[kk] = Buffer[kM] ^ (y >> 1) ^ ((y & 1u) * MATRIX_A);
            }
        }
        // 31-bit int [-0x80000000, 0x7fffffff]
        public int Int32() { return (int)this.UInt32(); }
        // 31-bit int [0, 0x7fffffff]
        public int Int31() { return (int)(this.UInt32() >> 1); }

        public uint UInt(uint value) { return this.UInt32() % value; }
        public int Int(int value) { return (int)UInt((uint)value); }

        public uint UIntExact(uint value)
        {
            while (true)
            {
                uint y = this.UInt32();
                if (y < (1LU << 32) / value * value) return y % value;
            }
        }
        public int IntExact(int value) { return (int)UIntExact((uint)value); }

        public double Double() { return Double32CO(); }
        // 32-bit double [0, 1]
        public double Double32CC() { return this.UInt32() * (1.0 / ((1L << 32) - 1)); }
        // 32-bit double [0, 1)
        public double Double32CO() { return this.UInt32() * (1.0 / (1L << 32)); }
        // 32-bit double (0, 1]
        public double Double32OC() { return (this.UInt32() + 1.0) * (1.0 / (1L << 32)); }
        // 32-bit double (0, 1)
        public double Double32OO() { return (this.UInt32() + 0.5) * (1.0 / (1L << 32)); }
        // 53-bit double [0, 1]
        public double Double53CC() { return ((this.UInt32() >> 5) * (double)(1 << 26) + (this.UInt32() >> 6)) * (1.0 / ((1L << 53) - 1)); }
        // 53-bit double [0, 1)
        public double Double53CO() { return ((this.UInt32() >> 5) * (double)(1 << 26) + (this.UInt32() >> 6)) * (1.0 / (1L << 53)); }
        // 53-bit double (0, 1]
        public double Double53OC() { return ((this.UInt32() >> 5) * (double)(1 << 26) + (this.UInt32() >> 6) + 1.0) * (1.0 / (1L << 53)); }
        // 52-bit double (0, 1)
        public double Double52OO() { return ((this.UInt32() >> 6) * (double)(1 << 26) + (this.UInt32() >> 6) + 0.5) * (1.0 / (1L << 52)); }

        public double Gaussian() { return Gaussian32(); }
        double BufferGaussian32 = double.PositiveInfinity;
        public double Gaussian32()
        {
            if (BufferGaussian32 != double.PositiveInfinity) { double x = BufferGaussian32; BufferGaussian32 = double.PositiveInfinity; return x; }
            double v1, v2, r;
            do
            {
                v1 = this.Double32OO() - 0.5;
                v2 = this.Double32OO() - 0.5;
                r = v1 * v1 + v2 * v2;
            } while (r >= 0.25 || r == 0);
            double f = Math.Sqrt(-2 * Math.Log(r * 4) / r);
            BufferGaussian32 = v2 * f;
            return v1 * f;
        }
        double BufferGaussian52 = double.PositiveInfinity;
        public double Gaussian52()
        {
            if (BufferGaussian52 != double.PositiveInfinity) { double x = BufferGaussian52; BufferGaussian52 = double.PositiveInfinity; return x; }
            double v1, v2, r;
            do
            {
                v1 = this.Double52OO() - 0.5;
                v2 = this.Double52OO() - 0.5;
                r = v1 * v1 + v2 * v2;
            } while (r >= 0.25 || r == 0);
            double f = Math.Sqrt(-2 * Math.Log(r * 4) / r);
            BufferGaussian52 = v2 * f;
            return v1 * f;
        }
        public int Categorical(double[] distribution)
        {
            double r = Double();
            double c = 0.0;
            for (int i = 0; i < distribution.Length - 1; i++)
            {
                c += distribution[i];
                if (r < c) return i;
            }
            return distribution.Length - 1;
        }
    }
    #endregion

    #region Comparer classes
    public class EqualityComparerArray<T> : IEqualityComparer<T[]>
        where T : IComparable<T>
    {
        public EqualityComparerArray() { }
        #region IEqualityComparer<T[]> メンバー
        public bool Equals(T[] x, T[] y) { return Ex.Compare(x, y) == 0; }
        public int GetHashCode(T[] obj) { int a = 0; foreach (var item in obj) a += item.GetHashCode() * 3; return a; }
        #endregion
    }
    public class EqualityComparerArray2<T> : IEqualityComparer<T[][]>
        where T : IComparable<T>
    {
        public EqualityComparerArray2() { }
        #region IEqualityComparer<T[][]> メンバー
        public bool Equals(T[][] x, T[][] y) { return Ex.Compare(x, y) == 0; }
        public int GetHashCode(T[][] obj) { int a = 0; foreach (var item in obj) a += item.GetHashCode() * 3; return a; }
        #endregion
    }
    public class EqualityComparerArray3<T> : IEqualityComparer<T[][][]>
        where T : IComparable<T>
    {
        public EqualityComparerArray3() { }
        #region IEqualityComparer<T[][][]> メンバー
        public bool Equals(T[][][] x, T[][][] y) { return Ex.Compare(x, y) == 0; }
        public int GetHashCode(T[][][] obj) { int a = 0; foreach (var item in obj) a += item.GetHashCode() * 3; return a; }
        #endregion
    }
    public class EqualityComparerArray4<T> : IEqualityComparer<T[][][][]>
        where T : IComparable<T>
    {
        public EqualityComparerArray4() { }
        #region IEqualityComparer<T[][][][]> メンバー
        public bool Equals(T[][][][] x, T[][][][] y) { return Ex.Compare(x, y) == 0; }
        public int GetHashCode(T[][][][] obj) { int a = 0; foreach (var item in obj) a += item.GetHashCode() * 3; return a; }
        #endregion
    }
    public class EqualityComparerIList<T> : IEqualityComparer<IList<T>>
        where T : IComparable<T>
    {
        public EqualityComparerIList() { }
        #region IEqualityComparer<IList<T>> メンバー
        public bool Equals(IList<T> x, IList<T> y) { return Ex.Compare(x, y) == 0; }
        public int GetHashCode(IList<T> obj) { int a = 0; foreach (var item in obj) a += item.GetHashCode() * 3; return a; }
        #endregion
    }

    public class ComparerArray<T> : IComparer<T[]>
        where T : IComparable<T>
    {
        public ComparerArray() { }
        #region IComparer<T[]> メンバー
        public int Compare(T[] x, T[] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerArray2<T> : IComparer<T[][]>
        where T : IComparable<T>
    {
        public ComparerArray2() { }
        #region IComparer<T[][]> メンバー
        public int Compare(T[][] x, T[][] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerArray3<T> : IComparer<T[][][]>
        where T : IComparable<T>
    {
        public ComparerArray3() { }
        #region IComparer<T[][][]> メンバー
        public int Compare(T[][][] x, T[][][] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerArray4<T> : IComparer<T[][][][]>
        where T : IComparable<T>
    {
        public ComparerArray4() { }
        #region IComparer<T[][][]> メンバー
        public int Compare(T[][][][] x, T[][][][] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerIList<T> : IComparer<IList<T>>
        where T : IComparable<T>
    {
        public ComparerIList() { }
        #region IComparer<IList<T>> メンバー
        public int Compare(IList<T> x, IList<T> y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerComparison<T> : IComparer<T>
    {
        Comparison<T> comparer;
        public ComparerComparison(Comparison<T> comparison) { this.comparer = comparison; }
        #region IComparer<T> メンバ
        public int Compare(T x, T y) { return comparer(x, y); }
        #endregion
    }
    public class ComparerReverse<T> : IComparer<T>
    {
        IComparer<T> comparer;
        public ComparerReverse() { this.comparer = Comparer<T>.Default; }
        public ComparerReverse(IComparer<T> comparer) { this.comparer = comparer; }
        public ComparerReverse(Comparison<T> comparer) { this.comparer = new ComparerComparison<T>(comparer); }
        #region IComparer<T> メンバ
        public int Compare(T x, T y) { return comparer.Compare(y, x); }
        #endregion
    }
    #endregion

    #region List classes
    [Serializable]
    public class SortedList<T> : List<T>
    {
        protected readonly IComparer<T> _Comparer;

        public SortedList()
            : base()
        {
            this._Comparer = Comparer<T>.Default;
        }
        public SortedList(IComparer<T> comparer)
            : base()
        {
            this._Comparer = comparer;
        }
        public SortedList(int capacity)
            : base(capacity)
        {
            this._Comparer = Comparer<T>.Default;
        }
        public SortedList(IComparer<T> comparer, int capacity)
            : base(capacity)
        {
            this._Comparer = comparer;
        }

        public IComparer<T> Comparer
        {
            get { return _Comparer; }
        }

        public bool AddOrDiscard(T item)
        {
            int index = IndexOf(item);
            if (index >= 0) return false;
            else { Insert(~index, item); return true; }
        }
        public bool AddOrOverwrite(T item)
        {
            int index = IndexOf(item);
            if (index >= 0) { this[index] = item; return false; }
            else { Insert(~index, item); return true; }
        }
        public T Pop()
        {
            T item = this[Count - 1];
            RemoveAt(Count - 1);
            return item;
        }

        public new int IndexOf(T item) { return BinarySearch(item); }
        public new void Add(T item)
        {
            int index = BinarySearch(item);
            if (index >= 0) { ThrowException.Argument("item"); return; }
            Insert(~index, item);
        }
        public new bool Contains(T item) { return BinarySearch(item) >= 0; }
    }

    [Serializable]
    public class ListedList<T> : IList<T>
    {
        protected const int FixedSize = 1024, FixedHalfSize = FixedSize / 2;
        protected List<List<T>> _Items = new List<List<T>>();
        protected int _Count = 0;

        public ListedList()
        {
        }
        public ListedList(IEnumerable<T> collections)
        {
            foreach (var item in collections) Add(item);
        }

        protected Int2 DecomposeIndex(int index)
        {
            int j = index;
            for (int i = 0; i < _Items.Count; i++)
            {
                if (j < _Items[i].Count) return new Int2(i, j);
                j -= _Items[i].Count;
            }
            return ThrowException.ArgumentOutOfRange<Int2>("index");
        }
        protected int ComposeIndex(Int2 indexes)
        {
            int index = indexes.Y;
            for (int i = indexes.X; --i >= 0; )
                index += _Items[i].Count;
            return index;
        }
        protected Int2 _IndexOf(T item)
        {
            for (int i = 0; i < _Items.Count; ++i)
            {
                List<T> items = _Items[i];
                int j = _Items[i].IndexOf(item);
                if (j >= 0) return new Int2(i, j);
            }
            return ~new Int2(0, 0);
        }
        protected Int2 _BinarySearch(T item) { return _BinarySearch(item, Comparer<T>.Default); }
        protected Int2 _BinarySearch(T item, IComparer<T> comparer)
        {
            int i0 = 0, i1 = _Items.Count;
            while (i0 < i1)
            {
                int i = (i0 + i1) / 2;
                int c = comparer.Compare(item, _Items[i][0]);
                if (c == 0) return new Int2(i, 0);
                if (c < 0) i1 = i; else i0 = i + 1;
            }
            if (i0 == 0) return ~new Int2(0, 0);
            int ii = i0 - 1;
            List<T> items = _Items[ii];
            int j0 = 0, j1 = items.Count;
            while (j0 < j1)
            {
                int j = (j0 + j1) / 2;
                int c = comparer.Compare(item, items[j]);
                if (c == 0) return new Int2(ii, j);
                if (c < 0) j1 = j; else j0 = j + 1;
            }
            return ~new Int2(ii, j0);
        }

        protected void _Insert(Int2 indexes, T item)
        {
            if (indexes.X == _Items.Count) _Items.Add(new List<T>(FixedSize));
            if (_Items[indexes.X].Count == FixedSize)
            {
                _Items.Insert(indexes.X + 1, new List<T>(FixedSize));
                _Items[indexes.X + 1].AddRange(_Items[indexes.X].GetRange(FixedHalfSize, FixedHalfSize));
                _Items[indexes.X].RemoveRange(FixedHalfSize, FixedHalfSize);
                if (indexes.Y >= FixedHalfSize) indexes -= new Int2(-1, FixedHalfSize);
            }
            _Items[indexes.X].Insert(indexes.Y, item);
            _Count++;
        }
        protected void _RemoveAt(Int2 indexes)
        {
            _Count--;
            _Items[indexes.X].RemoveAt(indexes.Y);
            if (_Items[indexes.X].Count == 0) _Items.RemoveAt(indexes.X);
        }

        public int RemoveAll(Predicate<T> match)
        {
            int count = 0;
            _Items.RemoveAll(items =>
            {
                count += items.RemoveAll(match);
                return items.Count == 0;
            });
            _Count -= count;
            return count;
        }
        public void CopyTo(T[] array)
        {
            int index = 0;
            foreach (List<T> items in _Items) { items.CopyTo(array, index); index += items.Count; }
        }
        public IEnumerable<T> Reverse()
        {
            for (int i = _Items.Count; --i >= 0; )
            {
                List<T> items = _Items[i];
                for (int j = items.Count; --j >= 0; )
                    yield return items[j];
            }
        }
        public override string ToString()
        {
            return string.Format("Count = {0}", _Count);
        }

        public T Pop()
        {
            _Count--;
            List<T> items = _Items[_Items.Count - 1];
            T item = items[items.Count - 1];
            items.RemoveAt(items.Count - 1);
            if (items.Count == 0) _Items.RemoveAt(_Items.Count - 1);
            return item;
        }

        public void Sort()
        {
            var newlist = new ListedList<T>();
            foreach (var list in _Items) { list.Sort(); list.Reverse(); }
            var head = _Items.Select(list => list.Count - 1).ToList();
            var index = Enumerable.Range(0, _Items.Count).ToList();
            var comparer = new ComparerComparison<int>((y, x) => Comparer<T>.Default.Compare(_Items[x][head[x]], _Items[y][head[y]]));
            index.Sort(comparer);
            while (_Items.Count > 0)
            {
                int most = index.Pop();
                newlist.Add(_Items[most][head[most]]);
                if (--head[most] >= 0)
                {
                    int insert = ~index.BinarySearch(most, comparer);
                    if (insert < 0) insert = ~insert;
                    index.Insert(insert, most);
                }
                else
                {
                    for (int i = index.Count; --i >= 0; )
                        if (index[i] > most) index[i]--;
                    head.RemoveAt(most);
                    _Items.RemoveAt(most);
                }
            }
            _Items = newlist._Items;
        }

        public int LetDistinct()
        {
            var comparer = EqualityComparer<T>.Default;
            int overlap = 0;
            for (int i = _Items.Count; --i >= 0; )
            {
                var list = _Items[i];
                for (int j = list.Count; --j > 0; )
                {
                    if (comparer.Equals(list[j], list[j - 1])) { list.RemoveAt(j); overlap++; }
                }
                if (i > 0)
                {
                    if (comparer.Equals(list[0], _Items[i - 1].Last())) { list.RemoveAt(0); overlap++; }
                }
                if (list.Count == 0) _Items.RemoveAt(i);
            }
            _Count -= overlap;
            return overlap;
        }

        #region IList<T> メンバ
        public int IndexOf(T item)
        {
            Int2 indices = _IndexOf(item);
            if (indices.X >= 0) return ComposeIndex(indices);
            else return ~ComposeIndex(~indices);
        }
        public void Insert(int index, T item) { _Insert(DecomposeIndex(index), item); }
        public void RemoveAt(int index) { _RemoveAt(DecomposeIndex(index)); }
        public T this[int index]
        {
            get { Int2 indices = DecomposeIndex(index); return _Items[indices.X][indices.Y]; }
            set { Int2 indices = DecomposeIndex(index); _Items[indices.X][indices.Y] = value; }
        }
        #endregion
        #region ICollection<T> メンバ
        public void Add(T item)
        {
            if (_Items.Count == 0 || _Items[_Items.Count - 1].Count == FixedSize) _Items.Add(new List<T>(FixedSize));
            _Items[_Items.Count - 1].Add(item);
            _Count++;
        }
        public void Clear()
        {
            _Count = 0;
            _Items.Clear();
        }
        public bool Contains(T item)
        {
            return _IndexOf(item).X >= 0;
        }
        public void CopyTo(T[] array, int arrayIndex)
        {
            int index = arrayIndex;
            foreach (var items in _Items)
            {
                items.CopyTo(array, index);
                index += items.Count;
            }
        }
        public int Count
        {
            get { return _Count; }
        }
        public bool IsReadOnly
        {
            get { return false; }
        }
        public bool Remove(T item)
        {
            Int2 indices = _IndexOf(item);
            if (indices.X >= 0) { _RemoveAt(indices); return true; }
            else { return false; }
        }
        #endregion
        #region IEnumerable<T> メンバ
        public IEnumerator<T> GetEnumerator()
        {
            for (int i = 0; i < _Items.Count; ++i)
            {
                List<T> items = _Items[i];
                for (int j = 0; j < items.Count; ++j)
                    yield return items[j];
            }
        }
        #endregion
        #region IEnumerable メンバ
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
        #endregion
    }

    [Serializable]
    public class SortedListedList<T> : ListedList<T>
    {
        protected readonly IComparer<T> _Comparer;

        public SortedListedList()
            : base()
        {
            this._Comparer = Comparer<T>.Default;
        }
        public SortedListedList(IEnumerable<T> collections)
        {
            foreach (var item in collections) Add(item);
        }
        public SortedListedList(IComparer<T> comparer)
            : base()
        {
            this._Comparer = comparer;
        }

        public IComparer<T> Comparer
        {
            get { return _Comparer; }
        }

        protected new Int2 _BinarySearch(T item) { return _BinarySearch(item, _Comparer); }

        public int BinarySearch(T item) { return ComposeIndex(_BinarySearch(item)); }
        public new bool Add(T item)
        {
            Int2 indices = _BinarySearch(item);
            if (indices.X >= 0) { return false; }
            else { _Insert(~indices, item); return true; }
        }
        public new void Insert(int index, T item) { ThrowException.NotImplemented(); }

        public bool AddOrDiscard(T item)
        {
            Int2 indices = _BinarySearch(item);
            if (indices.X >= 0) return false;
            else { _Insert(~indices, item); return true; }
        }
        public bool AddOrOverwrite(T item)
        {
            Int2 indices = _BinarySearch(item);
            if (indices.X >= 0) { _Items[indices.X][indices.Y] = item; return false; }
            else { _Insert(~indices, item); return true; }
        }
        public T FindOrDefault(T item)
        {
            Int2 index = _BinarySearch(item);
            return index.X >= 0 ? _Items[index.X][index.Y] : default(T);
        }

        public new bool Contains(T item)
        {
            return _BinarySearch(item).X >= 0;
        }
    }
    #endregion

    // Extended constructors
    public static partial class New
    {
        #region miscs
        public static T Let<T>(ref T obj) where T : new() { return obj = new T(); }
        public static T Type<T>() where T : new() { return new T(); }
        public static TResult Func<TResult>(Func<TResult> func) { return func(); }
        public static KeyValuePair<TKey, TValue> KeyValuePair<TKey, TValue>(TKey key, TValue value) { return new KeyValuePair<TKey, TValue>(key, value); }
        public static V2<T0, T1> V2<T0, T1>(T0 v0, T1 v1) { return new V2<T0, T1>(v0, v1); }
        public static V3<T0, T1, T2> V3<T0, T1, T2>(T0 v0, T1 v1, T2 v2) { return new V3<T0, T1, T2>(v0, v1, v2); }
        public static V4<T0, T1, T2, T3> V4<T0, T1, T2, T3>(T0 v0, T1 v1, T2 v2, T3 v3) { return new V4<T0, T1, T2, T3>(v0, v1, v2, v3); }
        public static Array1<T> Array1<T>(T v0) { return new Array1<T>(v0); }
        public static Array2<T> Array2<T>(T v0, T v1) { return new Array2<T>(v0, v1); }
        public static Array3<T> Array3<T>(T v0, T v1, T v2) { return new Array3<T>(v0, v1, v2); }
        public static Array4<T> Array4<T>(T v0, T v1, T v2, T v3) { return new Array4<T>(v0, v1, v2, v3); }

        public static ParallelOptions ParallelOptions(int degree = 0)
        {
            var result = new ParallelOptions();
            result.MaxDegreeOfParallelism = (degree > 0) ? degree : Environment.ProcessorCount + degree;
            return result;
        }
        #endregion

        #region Array
        public static T[] Array<T>(int length, Func<int, T> selector)
        {
            var result = new T[length];
            for (int i = 0; i < length; i++)
                result[i] = selector(i);
            return result;
        }
        public static T[,] Array<T>(int length0, int length1, Func<int, int, T> selector)
        {
            var result = new T[length0, length1];
            for (int i0 = 0; i0 < length0; i0++)
                for (int i1 = 0; i1 < length1; i1++)
                    result[i0, i1] = selector(i0, i1);
            return result;
        }
        public static T[, ,] Array<T>(int length0, int length1, int length2, Func<int, int, int, T> selector)
        {
            var result = new T[length0, length1, length2];
            for (int i0 = 0; i0 < length0; i0++)
                for (int i1 = 0; i1 < length1; i1++)
                    for (int i2 = 0; i2 < length2; i2++)
                        result[i0, i1, i2] = selector(i0, i1, i2);
            return result;
        }
        public static T[,] Array<T>(Int2 lengths, Func<int, int, T> selector) { return Array(lengths.v0, lengths.v1, selector); }
        public static T[, ,] Array<T>(Int3 lengths, Func<int, int, int, T> selector) { return Array(lengths.v0, lengths.v1, lengths.v2, selector); }
        public static T[] ParallelArray<T>(int length, Func<int, T> selector)
        {
            var result = new T[length];
            Parallel.For(0, length, i => result[i] = selector(i));
            return result;
        }
        public static T[,] ParallelArray<T>(int length0, int length1, Func<int, int, T> selector)
        {
            var result = new T[length0, length1];
            Parallel.ForEach(Ex.Range(length0, length1), i => result[i.v0, i.v1] = selector(i.v0, i.v1));
            return result;
        }
        public static T[, ,] ParallelArray<T>(int length0, int length1, int length2, Func<int, int, int, T> selector)
        {
            var result = new T[length0, length1, length2];
            Parallel.ForEach(Ex.Range(length0, length1, length2), i => result[i.v0, i.v1, i.v2] = selector(i.v0, i.v1, i.v2));
            return result;
        }

        public static List<T> List<T>(int length, Func<int, T> selector)
        {
            var result = new List<T>(length);
            for (int i = 0; i < length; i++)
                result.Add(selector(i));
            return result;
        }

        public static T[] Array<T>(int length, T item)
        {
            var result = new T[length];
            for (int i = 0; i < length; i++)
                result[i] = item;
            return result;
        }
        public static T[,] Array<T>(int length0, int length1, T item)
        {
            var result = new T[length0, length1];
            for (int i0 = 0; i0 < length0; i0++)
                for (int i1 = 0; i1 < length1; i1++)
                    result[i0, i1] = item;
            return result;
        }
        public static T[, ,] Array<T>(int length0, int length1, int length2, T item)
        {
            var result = new T[length0, length1, length2];
            for (int i0 = 0; i0 < length0; i0++)
                for (int i1 = 0; i1 < length1; i1++)
                    for (int i2 = 0; i2 < length2; i2++)
                        result[i0, i1, i2] = item;
            return result;
        }
        public static T[,] Array<T>(Int2 lengths, T item) { return Array(lengths.v0, lengths.v1, item); }
        public static T[, ,] Array<T>(Int3 lengths, T item) { return Array(lengths.v0, lengths.v1, lengths.v2, item); }
        public static List<T> List<T>(int length, T item)
        {
            var result = new List<T>(length);
            for (int i = 0; i < length; i++)
                result.Add(item);
            return result;
        }

        public static T[] Array<T>(int length) { return new T[length]; }
        public static T[,] Array<T>(int length0, int length1) { return new T[length0, length1]; }
        public static T[, ,] Array<T>(int length0, int length1, int length2) { return new T[length0, length1, length2]; }
        public static T[,] Array<T>(Int2 lengths) { return new T[lengths.v0, lengths.v1]; }
        public static T[, ,] Array<T>(Int3 lengths) { return new T[lengths.v0, lengths.v1, lengths.v2]; }
        #endregion

        #region Cache
        public static Func<T, TResult> Cache<T, TResult>(this Func<T, TResult> function)
        {
            var dictionary = new Dictionary<T, TResult>();
            return (arg) =>
            {
                lock (dictionary)
                {
                    if (!dictionary.ContainsKey(arg)) dictionary.Add(arg, function(arg));
                    return dictionary[arg];
                }
            };
        }
        public static Func<T1, T2, TResult> Cache<T1, T2, TResult>(this Func<T1, T2, TResult> function)
        {
            var dictionary = new Dictionary<V2<T1, T2>, TResult>();
            return (arg1, arg2) =>
            {
                var key = New.V2(arg1, arg2);
                lock (dictionary)
                {
                    if (!dictionary.ContainsKey(key)) dictionary.Add(key, function(arg1, arg2));
                    return dictionary[key];
                }
            };
        }
        public static Func<T1, T2, T3, TResult> Cache<T1, T2, T3, TResult>(this Func<T1, T2, T3, TResult> function)
        {
            var dictionary = new Dictionary<V3<T1, T2, T3>, TResult>();
            return (arg1, arg2, arg3) =>
            {
                var key = New.V3(arg1, arg2, arg3);
                lock (dictionary)
                {
                    if (!dictionary.ContainsKey(key)) dictionary.Add(key, function(arg1, arg2, arg3));
                    return dictionary[key];
                }
            };
        }
        public static Func<T1, T2, T3, T4, TResult> Cache<T1, T2, T3, T4, TResult>(this Func<T1, T2, T3, T4, TResult> function)
        {
            var dictionary = new Dictionary<V4<T1, T2, T3, T4>, TResult>();
            return (arg1, arg2, arg3, arg4) =>
            {
                var key = New.V4(arg1, arg2, arg3, arg4);
                lock (dictionary)
                {
                    if (!dictionary.ContainsKey(key)) dictionary.Add(key, function(arg1, arg2, arg3, arg4));
                    return dictionary[key];
                }
            };
        }
        #endregion

        #region Culture-based
        public static StreamReader StreamReader(string path)
        {
            var encoding = Ct.DetectEncoding(path);
            if (encoding == null) return null;
            return new StreamReader(path, encoding);
        }
        #endregion
    }

    // Utility and extension methods
    public static partial class Ex
    {
        #region Non-extension methods
        // System.Web.Util.HashCodeCombiner, System.Tuple
        public static int CombineHashCodes(int h1, int h2) { return (((h1 << 5) + h1) ^ h2); }
        public static int CombineHashCodes(int h1, int h2, int h3) { return CombineHashCodes(CombineHashCodes(h1, h2), h3); }
        public static int CombineHashCodes(int h1, int h2, int h3, int h4) { return CombineHashCodes(CombineHashCodes(h1, h2), CombineHashCodes(h3, h4)); }
        public static int CombineHashCodes(int h1, int h2, int h3, int h4, int h5) { return CombineHashCodes(CombineHashCodes(h1, h2, h3, h4), h5); }
        public static int CombineHashCodes(int h1, int h2, int h3, int h4, int h5, int h6) { return CombineHashCodes(CombineHashCodes(h1, h2, h3, h4), CombineHashCodes(h5, h6)); }
        public static int CombineHashCodes(int h1, int h2, int h3, int h4, int h5, int h6, int h7) { return CombineHashCodes(CombineHashCodes(h1, h2, h3, h4), CombineHashCodes(h5, h6, h7)); }
        public static int CombineHashCodes(int h1, int h2, int h3, int h4, int h5, int h6, int h7, int h8) { return CombineHashCodes(CombineHashCodes(h1, h2, h3, h4), CombineHashCodes(h5, h6, h7, h8)); }

        public static IEnumerable<int> Range(int count) { return Enumerable.Range(0, count); }
        public static IEnumerable<int> ParallelRange(int count) { return ParallelEnumerable.Range(0, count); }
        public static IEnumerable<Int2> Range(int count0, int count1)
        {
            var i = default(Int2);
            for (i.v0 = 0; i.v0 < count0; i.v0++)
                for (i.v1 = 0; i.v1 < count1; i.v1++)
                    yield return i;
        }
        public static IEnumerable<Int3> Range(int count0, int count1, int count2)
        {
            var i = default(Int3);
            for (i.v0 = 0; i.v0 < count0; i.v0++)
                for (i.v1 = 0; i.v1 < count1; i.v1++)
                    for (i.v2 = 0; i.v2 < count2; i.v2++)
                        yield return i;
        }
        public static IEnumerable<Int2> Range(Int2 counts) { return Range(counts.v0, counts.v1); }
        public static IEnumerable<Int3> Range(Int3 counts) { return Range(counts.v0, counts.v1, counts.v2); }

        public static IEnumerable<int[]> Range(this int[] counts)
        {
            int[] index = EnumeratorReset(counts);
            if (index == null) yield break;
            do
                yield return index;
            while (index.EnumeratorMoveNext(counts));
        }
        public static int[] EnumeratorReset(int[] counts)
        {
            if (counts == null || counts.Length == 0) return null;
            for (int i = 0; i < counts.Length; i++)
                if (counts[i] <= 0) return null;
            return new int[counts.Length];
        }
        public static bool EnumeratorMoveNext(this int[] index, int[] counts)
        {
            for (int i = index.Length; --i >= 0; )
            {
                if (++index[i] < counts[i]) return true;
                index[i] = 0;
            }
            return false;
        }

        public static IEnumerable<int> FromTo(int start, int end)
        {
            if (start <= end)
                for (int i = start; i < end; i++) yield return i;
            else
                for (int i = start; --i >= end; ) yield return i;
        }
        public static IEnumerable<int> FromToStep(int start, int end, int step)
        {
            if (start <= end)
                for (int i = start; i < end; i += step) yield return i;
            else
                for (int i = start; (i -= step) >= end; ) yield return i;
        }
        public static IEnumerable<int> FromUntil(int start, int end)
        {
            if (start <= end)
                for (int i = start; i <= end; i++) yield return i;
            else
                for (int i = start; i >= end; i--) yield return i;
        }
        public static IEnumerable<int> FromUntilStep(int start, int end, int step)
        {
            if (start <= end)
                for (int i = start; i <= end; i += step) yield return i;
            else
                for (int i = start; i >= end; i -= step) yield return i;
        }

        public static void Repeat(int count, Action action)
        {
            for (int i = 0; i < count; i++)
                action();
        }
        public static void For(int count, Action<int> action)
        {
            for (int i = 0; i < count; i++)
                action(i);
        }
        public static void For(int count0, int count1, Action<int, int> action)
        {
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    action(i0, i1);
        }
        public static void For(int count0, int count1, int count2, Action<int, int, int> action)
        {
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        action(i0, i1, i2);
        }
        public static void For(int count0, int count1, int count2, int count3, Action<int, int, int, int> action)
        {
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        for (int i3 = 0; i3 < count3; i3++)
                            action(i0, i1, i2, i3);
        }
        public static void For(Int2 counts, Action<int, int> action) { For(counts.v0, counts.v1, action); }
        public static void For(Int3 counts, Action<int, int, int> action) { For(counts.v0, counts.v1, counts.v2, action); }
        public static void ParallelFor(int count, Action<int> action)
        {
            Parallel.For(0, count, i => action(i));
        }
        public static void ParallelFor(int count0, int count1, Action<int, int> action)
        {
            Parallel.ForEach(Range(count0, count1), i => action(i.v0, i.v1));
        }
        public static void ParallelFor(int count0, int count1, int count2, Action<int, int, int> action)
        {
            Parallel.ForEach(Range(count0, count1, count2), i => action(i.v0, i.v1, i.v2));
        }

        public static void ForEach<TSource>(this IEnumerable<TSource> source, Action<TSource> action)
        {
            foreach (var element in source) action(element);
        }
        public static void ForEach<TSource>(this IEnumerable<TSource> source, Action<TSource, int> action)
        {
            int c = 0;
            foreach (var element in source) action(element, c++);
        }
        public static void ParallelForEach<TSource>(this IEnumerable<TSource> source, Action<TSource> action)
        {
            Parallel.ForEach(source, action);
        }

        public static byte[] SerializeToByteArray(object obj)
        {
            if (obj == null) return new byte[0];
            var stream = new MemoryStream();
            new BinaryFormatter().Serialize(stream, obj);
            return stream.ToArray();
        }
        public static object DeserializeFromByteArray(byte[] array)
        {
            if (array.Length == 0) return null;
            var stream = new MemoryStream(array);
            return new BinaryFormatter().Deserialize(stream);
        }
        #endregion

        #region IEnumerable
        public static TSource FirstOrDefault<TSource>(this IEnumerable<TSource> source, TSource defaultvalue)
        {
            foreach (var element in source)
                return element;
            return defaultvalue;
        }
        public static TSource FirstOrDefault<TSource>(this IEnumerable<TSource> source, Func<TSource, bool> predicate, TSource defaultvalue)
        {
            foreach (var element in source)
                if (predicate(element)) return element;
            return defaultvalue;
        }

        public static IOrderedEnumerable<TSource> OrderBy<TSource>(this IEnumerable<TSource> source) { return source.OrderBy(x => x); }
        public static IOrderedEnumerable<TSource> OrderByDescending<TSource>(this IEnumerable<TSource> source) { return source.OrderByDescending(x => x); }

        public static Dictionary<TKey, TValue> ToDictionary<TKey, TValue>(this IEnumerable<KeyValuePair<TKey, TValue>> source)
        {
            var dictionary = new Dictionary<TKey, TValue>();
            foreach (var pair in source)
                dictionary.Add(pair.Key, pair.Value);
            return dictionary;
        }

        public static Dictionary<TKey, int> ToHistogram<TKey>(this IEnumerable<TKey> source, IEqualityComparer<TKey> comparer = null)
        {
            var dictionary = comparer == null ? new Dictionary<TKey, int>() : new Dictionary<TKey, int>(comparer);
            foreach (var e in source)
                if (dictionary.ContainsKey(e)) dictionary[e]++;
                else dictionary.Add(e, 1);
            return dictionary;
        }

        #region MinMax
        public static int MinIndex<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MinIndexItem().v0; }
        public static int MaxIndex<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MaxIndexItem().v0; }
        public static Array2<int> MinMaxIndex<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { var r = source.MinMaxIndexItem(); return New.Array2(r.v0.v0, r.v1.v0); }
        public static TSource MinItem<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MinIndexItem().v1; }
        public static TSource MaxItem<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MaxIndexItem().v1; }
        public static Array2<TSource> MinMaxItem<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { var r = source.MinMaxIndexItem(); return New.Array2(r.v0.v1, r.v1.v1); }
        public static V2<int, TSource> MinIndexItem<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MinMaxIndexItem_(true, false).v0; }
        public static V2<int, TSource> MaxIndexItem<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MinMaxIndexItem_(false, true).v1; }
        public static Array2<V2<int, TSource>> MinMaxIndexItem<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MinMaxIndexItem_(true, true); }
        static Array2<V2<int, TSource>> MinMaxIndexItem_<TSource>(this IEnumerable<TSource> source, bool minSw, bool maxSw)
            where TSource : IComparable<TSource>
        {
            var minIndex = default(int);
            var maxIndex = default(int);
            var minItem = default(TSource);
            var maxItem = default(TSource);
            int index = 0;
            foreach (var item in source)
            {
                if (minSw && (index == 0 || minItem.CompareTo(item) > 0)) { minIndex = index; minItem = item; }
                if (maxSw && (index == 0 || maxItem.CompareTo(item) < 0)) { maxIndex = index; maxItem = item; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.Array2(New.V2(minIndex, minItem), New.V2(maxIndex, maxItem));
        }

        public static int MinIndex<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinIndexItemValue(selector).v0; }
        public static int MaxIndex<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MaxIndexItemValue(selector).v0; }
        public static Array2<int> MinMaxIndex<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { var r = source.MinMaxIndexItemValue(selector); return New.Array2(r.v0.v0, r.v1.v0); }
        public static TSource MinItem<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinIndexItemValue(selector).v1; }
        public static TSource MaxItem<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MaxIndexItemValue(selector).v1; }
        public static Array2<TSource> MinMaxItem<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { var r = source.MinMaxIndexItemValue(selector); return New.Array2(r.v0.v1, r.v1.v1); }
        public static TValue MinValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinIndexItemValue(selector).v2; }
        public static TValue MaxValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MaxIndexItemValue(selector).v2; }
        public static Array2<TValue> MinMaxValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { var r = source.MinMaxIndexItemValue(selector); return New.Array2(r.v0.v2, r.v1.v2); }
        public static V3<int, TSource, TValue> MinIndexItemValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinMaxIndexItemValue_(true, false, selector).v0; }
        public static V3<int, TSource, TValue> MaxIndexItemValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinMaxIndexItemValue_(false, true, selector).v1; }
        public static Array2<V3<int, TSource, TValue>> MinMaxIndexItemValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinMaxIndexItemValue_(true, true, selector); }
        static Array2<V3<int, TSource, TValue>> MinMaxIndexItemValue_<TSource, TValue>(this IEnumerable<TSource> source, bool minSw, bool maxSw, Func<TSource, TValue> selector)
            where TValue : IComparable<TValue>
        {
            var minIndex = default(int);
            var maxIndex = default(int);
            var minItem = default(TSource);
            var maxItem = default(TSource);
            var minValue = default(TValue);
            var maxValue = default(TValue);
            int index = 0;
            foreach (var item in source)
            {
                var value = selector(item);
                if (minSw && (index == 0 || minValue.CompareTo(value) > 0)) { minIndex = index; minItem = item; minValue = value; }
                if (maxSw && (index == 0 || maxValue.CompareTo(value) < 0)) { maxIndex = index; maxItem = item; maxValue = value; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.Array2(New.V3(minIndex, minItem, minValue), New.V3(maxIndex, maxItem, maxValue));
        }
        #endregion
        #endregion

        #region Array
        public static void SizeCheck(Array array0, Array array1)
        {
            if (array0.Rank != array1.Rank || Ex.Range(array0.Rank).Any(d => array0.GetLength(d) != array1.GetLength(d))) ThrowException.SizeMismatch();
        }
        public static void SizeCheck<TSource0, TSource1>(TSource0[] array0, TSource1[] array1)
        {
            if (array0.Length != array1.Length) ThrowException.SizeMismatch();
        }
        public static void SizeCheck<TSource0, TSource1>(TSource0[,] array0, TSource1[,] array1)
        {
            if (array0.GetLength(0) != array1.GetLength(0) ||
                array0.GetLength(1) != array1.GetLength(1)) ThrowException.SizeMismatch();
        }
        public static void SizeCheck<TSource0, TSource1>(TSource0[, ,] array0, TSource1[, ,] array1)
        {
            if (array0.GetLength(0) != array1.GetLength(0) ||
                array0.GetLength(1) != array1.GetLength(1) ||
                array0.GetLength(2) != array1.GetLength(2)) ThrowException.SizeMismatch();
        }
        public static void SizeCheck(Array array0, Array array1, Array array2) { SizeCheck(array0, array1); SizeCheck(array0, array2); }
        public static void SizeCheck<TSource0, TSource1, TSource2>(TSource0[] array0, TSource1[] array1, TSource2[] array2) { SizeCheck(array0, array1); SizeCheck(array0, array2); }
        public static void SizeCheck<TSource0, TSource1, TSource2>(TSource0[,] array0, TSource1[,] array1, TSource2[,] array2) { SizeCheck(array0, array1); SizeCheck(array0, array2); }
        public static void SizeCheck<TSource0, TSource1, TSource2>(TSource0[, ,] array0, TSource1[, ,] array1, TSource2[, ,] array2) { SizeCheck(array0, array1); SizeCheck(array0, array2); }

        public static int IndexOf<T>(this T[] array, T value) { return Array.IndexOf(array, value); }
        public static void Clear(this Array array) { Array.Clear(array, 0, array.Length); }
        public static void CopyTo(this Array sourceArray, Array destinationArray) { sourceArray.CopyTo(destinationArray, 0); }
        public static void CopyTo(this Array sourceArray, int srcIndex, Array destinationArray, int dstinationIndex, int length) { Array.Copy(sourceArray, srcIndex, destinationArray, dstinationIndex, length); }
        public static T[] CloneX<T>(this T[] array) { return (T[])array.Clone(); }
        public static T[,] CloneX<T>(this T[,] array) { return (T[,])array.Clone(); }
        public static T[, ,] CloneX<T>(this T[, ,] array) { return (T[, ,])array.Clone(); }
        public static T[][] DeepClone<T>(this T[][] array) { return New.Array(array.Length, i => array[i].CloneX()); }
        public static T[][][] DeepClone<T>(this T[][][] array) { return New.Array(array.Length, i => array[i].DeepClone()); }

        public static int Lengths<T>(this T[] array) { return array.Length; }
        public static Int2 Lengths<T>(this T[,] array) { return new Int2(array.GetLength(0), array.GetLength(1)); }
        public static Int3 Lengths<T>(this T[, ,] array) { return new Int3(array.GetLength(0), array.GetLength(1), array.GetLength(2)); }
        public static int[] GetLengths<T>(this T[] array) { return new[] { array.Length }; }
        public static int[] GetLengths<T>(this T[,] array) { return new[] { array.GetLength(0), array.GetLength(1) }; }
        public static int[] GetLengths<T>(this T[, ,] array) { return new[] { array.GetLength(0), array.GetLength(1), array.GetLength(2) }; }
        public static int[] GetLengths(this Array array) { return New.Array(array.Rank, i => array.GetLength(i)); }

        public static IEnumerable<T> ToEnumerable<T>(this T[,] array)
        {
            var l = array.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    yield return array[i0, i1];
        }
        public static IEnumerable<T> ToEnumerable<T>(this T[, ,] array)
        {
            var l = array.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    for (int i2 = 0; i2 < l.v2; i2++)
                        yield return array[i0, i1, i2];
        }
        public static T[] ToArray<T>(this T[,] array)
        {
            var a = new T[array.Length];
            int c = 0;
            foreach (var item in array) a[c++] = item;
            return a;
        }
        public static T[] ToArray<T>(this T[, ,] array)
        {
            var a = new T[array.Length];
            int c = 0;
            foreach (var item in array) a[c++] = item;
            return a;
        }
        public static T[,] ToArray<T>(this IEnumerable<T> source, Int2 lengths)
        {
            var array = new T[lengths.v0, lengths.v1];
            int i0 = 0, i1 = 0;
            foreach (var element in source)
            {
                array[i0, i1] = element;
                if (++i1 == lengths.v1) { i1 = 0; ++i0; }
            }
            if (i0 != lengths.v0 || i1 != 0) ThrowException.InvalidOperation("not enough items");
            return array;
        }
        public static Array ToArray<T>(this IEnumerable<T> source, int[] lengths)
        {
            var index = EnumeratorReset(lengths);
            if (index == null) ThrowException.ArgumentOutOfRange("lengths");
            Array array = Array.CreateInstance(typeof(T), lengths);
            foreach (var element in source)
            {
                array.SetValue(element, index);
                if (!index.EnumeratorMoveNext(lengths)) break;
            }
            return array;
        }

        public static T[] Let<T>(this T[] array, Func<T, T> selector)
        {
            var l = array.Length;
            for (int i = 0; i < l; i++)
                array[i] = selector(array[i]);
            return array;
        }
        public static T[,] Let<T>(this T[,] array, Func<T, T> selector)
        {
            var l = array.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    array[i0, i1] = selector(array[i0, i1]);
            return array;
        }
        public static T[, ,] Let<T>(this T[, ,] array, Func<T, T> selector)
        {
            var l = array.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    for (int i2 = 0; i2 < l.v2; i2++)
                        array[i0, i1, i2] = selector(array[i0, i1, i2]);
            return array;
        }

        public static T[] Let<T>(this T[] array, Func<T, int, T> selector)
        {
            var l = array.Length;
            for (int i = 0; i < l; i++)
                array[i] = selector(array[i], i);
            return array;
        }
        public static T[,] Let<T>(this T[,] array, Func<T, int, int, T> selector)
        {
            var l = array.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    array[i0, i1] = selector(array[i0, i1], i0, i1);
            return array;
        }
        public static T[, ,] Let<T>(this T[, ,] array, Func<T, int, int, int, T> selector)
        {
            var l = array.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    for (int i2 = 0; i2 < l.v2; i2++)
                        array[i0, i1, i2] = selector(array[i0, i1, i2], i0, i1, i2);
            return array;
        }

        public static T[][] ToJaggedArray<T>(this T[,] array)
        {
            var result = new T[array.GetLength(0)][];
            for (int i = 0; i < result.Length; i++)
            {
                var row = new T[array.GetLength(1)];
                for (int j = 0; j < row.Length; j++)
                    row[j] = array[i, j];
                result[i] = row;
            }
            return result;
        }

        public static T[] ZeroTo<T>(this T[] array) { return new T[array.Length]; }
        public static T[,] ZeroTo<T>(this T[,] array) { return new T[array.GetLength(0), array.GetLength(1)]; }
        public static T[, ,] ZeroTo<T>(this T[, ,] array) { return new T[array.GetLength(0), array.GetLength(1), array.GetLength(2)]; }
        public static TResult[] ZeroTo<TSource, TResult>(this TSource[] array) { return new TResult[array.Length]; }
        public static TResult[,] ZeroTo<TSource, TResult>(this TSource[,] array) { return new TResult[array.GetLength(0), array.GetLength(1)]; }
        public static TResult[, ,] ZeroTo<TSource, TResult>(this TSource[, ,] array) { return new TResult[array.GetLength(0), array.GetLength(1), array.GetLength(2)]; }

        public static TResult[] SelectTo<T, TResult>(this IEnumerable<T> source, Func<T, TResult> selector) { return source.Select(selector).ToArray(); }
        public static TResult[] SelectTo<T, TResult>(this T[] array, Func<T, TResult> selector) { return New.Array(array.Length, i => selector(array[i])); }
        public static TResult[,] SelectTo<T, TResult>(this T[,] array, Func<T, TResult> selector) { return New.Array(array.GetLength(0), array.GetLength(1), (i0, i1) => selector(array[i0, i1])); }
        public static TResult[, ,] SelectTo<T, TResult>(this T[, ,] array, Func<T, TResult> selector) { return New.Array(array.GetLength(0), array.GetLength(1), array.GetLength(2), (i0, i1, i2) => selector(array[i0, i1, i2])); }
        public static TResult[] ParallelSelectTo<T, TResult>(this T[] array, Func<T, TResult> selector) { return New.ParallelArray(array.Length, i => selector(array[i])); }
        public static TResult[,] ParallelSelectTo<T, TResult>(this T[,] array, Func<T, TResult> selector) { return New.ParallelArray(array.GetLength(0), array.GetLength(1), (i0, i1) => selector(array[i0, i1])); }
        public static TResult[, ,] ParallelSelectTo<T, TResult>(this T[, ,] array, Func<T, TResult> selector) { return New.ParallelArray(array.GetLength(0), array.GetLength(1), array.GetLength(2), (i0, i1, i2) => selector(array[i0, i1, i2])); }

        public static TResult[] ZipTo<TSource0, TSource1, TResult>(this TSource0[] array0, TSource1[] array1, Func<TSource0, TSource1, TResult> selector) { Ex.SizeCheck(array0, array1); return New.Array(array0.Length, i => selector(array0[i], array1[i])); }
        public static TResult[,] ZipTo<TSource0, TSource1, TResult>(this TSource0[,] array0, TSource1[,] array1, Func<TSource0, TSource1, TResult> selector) { Ex.SizeCheck(array0, array1); return New.Array(array0.GetLength(0), array0.GetLength(1), (i0, i1) => selector(array0[i0, i1], array1[i0, i1])); }
        public static TResult[, ,] ZipTo<TSource0, TSource1, TResult>(this TSource0[, ,] array0, TSource1[, ,] array1, Func<TSource0, TSource1, TResult> selector) { Ex.SizeCheck(array0, array1); return New.Array(array0.GetLength(0), array0.GetLength(1), array0.GetLength(2), (i0, i1, i2) => selector(array0[i0, i1, i2], array1[i0, i1, i2])); }
        public static TResult[] ParallelZipTo<TSource0, TSource1, TResult>(this TSource0[] array0, TSource1[] array1, Func<TSource0, TSource1, TResult> selector) { Ex.SizeCheck(array0, array1); return New.ParallelArray(array0.Length, i => selector(array0[i], array1[i])); }
        public static TResult[,] ParallelZipTo<TSource0, TSource1, TResult>(this TSource0[,] array0, TSource1[,] array1, Func<TSource0, TSource1, TResult> selector) { Ex.SizeCheck(array0, array1); return New.ParallelArray(array0.GetLength(0), array0.GetLength(1), (i0, i1) => selector(array0[i0, i1], array1[i0, i1])); }
        public static TResult[, ,] ParallelZipTo<TSource0, TSource1, TResult>(this TSource0[, ,] array0, TSource1[, ,] array1, Func<TSource0, TSource1, TResult> selector) { Ex.SizeCheck(array0, array1); return New.ParallelArray(array0.GetLength(0), array0.GetLength(1), array0.GetLength(2), (i0, i1, i2) => selector(array0[i0, i1, i2], array1[i0, i1, i2])); }

        public static TResult[] ZipTo<TSource0, TSource1, TSource2, TResult>(this TSource0[] array0, TSource1[] array1, TSource2[] array2, Func<TSource0, TSource1, TSource2, TResult> selector) { Ex.SizeCheck(array0, array1, array2); return New.Array(array0.Length, i => selector(array0[i], array1[i], array2[i])); }
        public static TResult[,] ZipTo<TSource0, TSource1, TSource2, TResult>(this TSource0[,] array0, TSource1[,] array1, TSource2[,] array2, Func<TSource0, TSource1, TSource2, TResult> selector) { Ex.SizeCheck(array0, array1, array2); return New.Array(array0.GetLength(0), array0.GetLength(1), (i0, i1) => selector(array0[i0, i1], array1[i0, i1], array2[i0, i1])); }
        public static TResult[, ,] ZipTo<TSource0, TSource1, TSource2, TResult>(this TSource0[, ,] array0, TSource1[, ,] array1, TSource2[, ,] array2, Func<TSource0, TSource1, TSource2, TResult> selector) { Ex.SizeCheck(array0, array1, array2); return New.Array(array0.GetLength(0), array0.GetLength(1), array0.GetLength(2), (i0, i1, i2) => selector(array0[i0, i1, i2], array1[i0, i1, i2], array2[i0, i1, i2])); }
        public static TResult[] ParallelZipTo<TSource0, TSource1, TSource2, TResult>(this TSource0[] array0, TSource1[] array1, TSource2[] array2, Func<TSource0, TSource1, TSource2, TResult> selector) { Ex.SizeCheck(array0, array1, array2); return New.ParallelArray(array0.Length, i => selector(array0[i], array1[i], array2[i])); }
        public static TResult[,] ParallelZipTo<TSource0, TSource1, TSource2, TResult>(this TSource0[,] array0, TSource1[,] array1, TSource2[,] array2, Func<TSource0, TSource1, TSource2, TResult> selector) { Ex.SizeCheck(array0, array1, array2); return New.ParallelArray(array0.GetLength(0), array0.GetLength(1), (i0, i1) => selector(array0[i0, i1], array1[i0, i1], array2[i0, i1])); }
        public static TResult[, ,] ParallelZipTo<TSource0, TSource1, TSource2, TResult>(this TSource0[, ,] array0, TSource1[, ,] array1, TSource2[, ,] array2, Func<TSource0, TSource1, TSource2, TResult> selector) { Ex.SizeCheck(array0, array1, array2); return New.ParallelArray(array0.GetLength(0), array0.GetLength(1), array0.GetLength(2), (i0, i1, i2) => selector(array0[i0, i1, i2], array1[i0, i1, i2], array2[i0, i1, i2])); }

        public static IEnumerable<V2<TSource0, TSource1>> Zip<TSource0, TSource1>(this TSource0[] array0, TSource1[] array1)
        {
            Ex.SizeCheck(array0, array1);
            var l = array0.Length;
            for (int i = 0; i < l; i++)
                yield return New.V2(array0[i], array1[i]);
        }
        public static IEnumerable<V2<TSource0, TSource1>> Zip<TSource0, TSource1>(this TSource0[,] array0, TSource1[,] array1)
        {
            Ex.SizeCheck(array0, array1);
            var l = array0.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    yield return New.V2(array0[i0, i1], array1[i0, i1]);
        }
        public static IEnumerable<V2<TSource0, TSource1>> Zip<TSource0, TSource1>(this TSource0[, ,] array0, TSource1[, ,] array1)
        {
            Ex.SizeCheck(array0, array1);
            var l = array0.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    for (int i2 = 0; i2 < l.v2; i2++)
                        yield return New.V2(array0[i0, i1, i2], array1[i0, i1, i2]);
        }
        public static IEnumerable<V3<TSource0, TSource1, TSource2>> Zip<TSource0, TSource1, TSource2>(this TSource0[] array0, TSource1[] array1, TSource2[] array2)
        {
            Ex.SizeCheck(array0, array1, array2);
            var l = array0.Length;
            for (int i = 0; i < l; i++)
                yield return New.V3(array0[i], array1[i], array2[i]);
        }
        public static IEnumerable<V3<TSource0, TSource1, TSource2>> Zip<TSource0, TSource1, TSource2>(this TSource0[,] array0, TSource1[,] array1, TSource2[,] array2)
        {
            Ex.SizeCheck(array0, array1, array2);
            var l = array0.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    yield return New.V3(array0[i0, i1], array1[i0, i1], array2[i0, i1]);
        }
        public static IEnumerable<V3<TSource0, TSource1, TSource2>> Zip<TSource0, TSource1, TSource2>(this TSource0[, ,] array0, TSource1[, ,] array1, TSource2[, ,] array2)
        {
            Ex.SizeCheck(array0, array1, array2);
            var l = array0.Lengths();
            for (int i0 = 0; i0 < l.v0; i0++)
                for (int i1 = 0; i1 < l.v1; i1++)
                    for (int i2 = 0; i2 < l.v2; i2++)
                        yield return New.V3(array0[i0, i1, i2], array1[i0, i1, i2], array2[i0, i1, i2]);
        }

        public static T[] SkipTo<T>(this T[] array, int start)
        {
            var result = new T[array.Length - start];
            Array.Copy(array, start, result, 0, result.Length);
            return result;
        }
        public static T[] SkipTakeTo<T>(this T[] array, int start, int count)
        {
            var result = new T[count];
            Array.Copy(array, start, result, 0, count);
            return result;
        }
        public static IEnumerable<T> Row<T>(this T[,] matrix, int row)
        {
            return Enumerable.Range(0, matrix.GetLength(1)).Select(j => matrix[row, j]);
        }
        public static IEnumerable<T> Col<T>(this T[,] matrix, int col)
        {
            return Enumerable.Range(0, matrix.GetLength(0)).Select(i => matrix[i, col]);
        }

        public static T[] ConcatTo<T>(this IEnumerable<T> source0, IEnumerable<T> source1) { return source0.Concat(source1).ToArray(); }
        public static T[] ConcatTo<T>(this IEnumerable<IEnumerable<T>> blocks)
        {
            var result = new List<T>();
            foreach (var b in blocks) result.AddRange(b);
            return result.ToArray();
        }
        public static T[,] ConcatTo<T>(this T[,][,] blocks)
        {
            var l0 = new int[blocks.GetLength(0)];
            var l1 = new int[blocks.GetLength(1)];
            for (int i0 = 0; i0 < l0.Length; i0++)
            {
                for (int i1 = 0; i1 < l1.Length; i1++)
                {
                    var b = blocks[i0, i1];
                    if (b == null) continue;
                    if (l0[i0] != 0 && l0[i0] != b.GetLength(0) || l1[i1] != 0 && l1[i1] != b.GetLength(1)) ThrowException.SizeMismatch();
                    l0[i0] = b.GetLength(0);
                    l1[i1] = b.GetLength(1);
                }
            }
            var r = new T[l0.Sum(), l1.Sum()];
            for (int o0 = 0, i0 = 0; i0 < l0.Length; o0 += l0[i0], i0++)
            {
                for (int o1 = 0, i1 = 0; i1 < l1.Length; o1 += l1[i1], i1++)
                {
                    var b = blocks[i0, i1];
                    if (b == null) continue;
                    Ex.For(l0[i0], j0 => Ex.For(l1[i1], j1 => r[j0 + o0, j1 + o1] = b[j0, j1]));
                }
            }
            return r;
        }
        #endregion

        #region List
        public static bool AddOrDiscard<T>(this List<T> list, T item)
        {
            if (list.Contains(item)) return false;
            list.Add(item); return true;
        }
        public static bool AddOrOverwrite<T>(this List<T> list, T item)
        {
            int index = list.IndexOf(item);
            if (index != -1) { list[index] = item; return false; }
            list.Add(item); return true;
        }
        #endregion

        #region IList
        public static int[] FindAllIndex<T>(this IList<T> list, Predicate<T> match)
        {
            var result = new List<int>();
            for (int i = 0; i < list.Count; i++)
                if (match(list[i])) result.Add(i);
            return result.ToArray();
        }
        public static T[] SkipTo<T>(this IList<T> list, int start)
        {
            var result = new T[list.Count - start];
            for (int i = 0; i < result.Length; i++)
                result[i] = list[i + start];
            return result;
        }
        public static T[] SkipTakeTo<T>(this IList<T> list, int start, int count)
        {
            var result = new T[count];
            for (int i = 0; i < count; i++)
                result[i] = list[i + start];
            return result;
        }
        public static bool AllBackward<T>(this IList<T> list, Predicate<T> match)
        {
            for (int i = list.Count; --i >= 0; )
                if (!match(list[i])) return false;
            return true;
        }
        public static bool AnyBackward<T>(this IList<T> list, Predicate<T> match)
        {
            for (int i = list.Count; --i >= 0; )
                if (match(list[i])) return false;
            return true;
        }

        public static int BinarySearch<T>(this IList<T> list, T item) where T : IComparable<T>
        {
            int i0 = 0;
            int i1 = list.Count;
            while (true)
            {
                if (i0 == i1) return ~i0;
                int ii = (i0 + i1) / 2;
                int c = item.CompareTo(list[ii]);
                if (c == 0) return ii;
                if (c < 0) i1 = ii; else i0 = ii + 1;
            }
        }
        public static int BinarySearch<T>(this IList<T> list, T item, IComparer<T> comparer)
        {
            int i0 = 0;
            int i1 = list.Count;
            while (true)
            {
                if (i0 == i1) return ~i0;
                int ii = (i0 + i1) / 2;
                int c = comparer.Compare(item, list[ii]);
                if (c == 0) return ii;
                if (c < 0) i1 = ii; else i0 = ii + 1;
            }
        }

        public static T First<T>(this IList<T> list) { return list[0]; }
        public static T Last<T>(this IList<T> list) { return list[list.Count - 1]; }
        public static T Pop<T>(this IList<T> list)
        {
            var item = list[list.Count - 1];
            list.RemoveAt(list.Count - 1);
            return item;
        }
        public static T Pull<T>(this IList<T> list, int index)
        {
            var item = list[index];
            list.RemoveAt(index);
            return item;
        }

        public static int LetDistinctSorted<T>(this IList<T> list)
        {
            var comparer = EqualityComparer<T>.Default;
            int i = 0;
            for (int j = 1; j < list.Count; j++)
            {
                if (!comparer.Equals(list[i], list[j]))
                    list[++i] = list[j];
            }
            i++;
            int overlap = list.Count - i;
            for (int j = list.Count; --j >= i; )
                list.RemoveAt(j);
            return overlap;
        }

        public static void Let<TSource>(this IList<TSource> list, Func<TSource, TSource> selector)
        {
            for (int i = 0; i < list.Count; i++)
                list[i] = selector(list[i]);
        }
        public static void Let<TSource>(this IList<TSource> list, Func<TSource, int, TSource> selector)
        {
            for (int i = 0; i < list.Count; i++)
                list[i] = selector(list[i], i);
        }
        #endregion

        #region ICollection, SortedSet, IDictionary, ISet
        // ICollection
        public static void AddNew<T>(this ICollection<T> source) where T : new()
        {
            source.Add(new T());
        }
        public static void AddNew<T>(this ICollection<T> source, int count) where T : new()
        {
            for (int i = count; --i >= 0; )
                source.Add(new T());
        }

        // LinkedList
        public static void Push<T>(this LinkedList<T> list, T value) { list.AddLast(value); }
        public static T Pop<T>(this LinkedList<T> list) { var r = list.Last.Value; list.RemoveLast(); return r; }
        public static void Enqueue<T>(this LinkedList<T> list, T value) { list.AddLast(value); }
        public static T Dequeue<T>(this LinkedList<T> list) { var r = list.First.Value; list.RemoveFirst(); return r; }

        // SortedSet
        public static T First<T>(this SortedSet<T> set) { return set.Min; }
        public static T Last<T>(this SortedSet<T> set) { return set.Max; }
        public static T Pop<T>(this SortedSet<T> set) { var r = set.Max; set.Remove(r); return r; }

        // IDictionary
        public static KeyValuePair<TKey, TValue> Pop<TKey, TValue>(this IDictionary<TKey, TValue> dictionary) { var r = dictionary.Last(); dictionary.Remove(r.Key); return r; }
        public static bool AddOrDiscard<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue value)
        {
            if (dictionary.ContainsKey(key)) return false;
            dictionary.Add(key, value); return true;
        }
        public static bool AddOrOverwrite<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue value)
        {
            if (dictionary.ContainsKey(key)) { dictionary[key] = value; return false; }
            dictionary.Add(key, value); return true;
        }
        public static TValue GetItemOrDefault<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key)
        {
            return dictionary.Keys.Contains(key) ? dictionary[key] : default(TValue);
        }
        public static TValue TryGetValue<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue defaultvalue)
        {
            TValue value;
            return dictionary.TryGetValue(key, out value) ? value : defaultvalue;
        }
        public static TValue TryGetValue<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key)
        {
            TValue value;
            return dictionary.TryGetValue(key, out value) ? value : default(TValue);
        }

        // ISet
        public static bool AddOrDiscard<T>(this ISet<T> set, T item)
        {
            return set.Add(item);
        }
        public static bool AddOrOverwrite<T>(this ISet<T> set, T item)
        {
            if (set.Contains(item)) { set.Remove(item); set.Add(item); return false; }
            else { set.Add(item); return true; }
        }

        // Dictionary
        // SortedDictionary
        // SortedList
        #endregion

        #region Sorting
        public static T[] LetReverse<T>(this T[] array) { Array.Reverse(array); return array; }
        public static T[] LetSort<T>(this T[] array) where T : IComparable<T> { Array.Sort(array); return array; }
        public static T[] LetSort<T>(this T[] array, Comparison<T> compare) { Array.Sort(array, compare); return array; }
        public static T[] LetSort<T>(this T[] array, Func<T, int> selector) { return array.LetSort((x, y) => selector(x).CompareTo(selector(y))); }
        public static T[] LetSort<T>(this T[] array, Func<T, double> selector) { return array.LetSort((x, y) => selector(x).CompareTo(selector(y))); }
        public static T[] LetSort<T>(this T[] array, Func<T, BigInteger> selector) { return array.LetSort((x, y) => selector(x).CompareTo(selector(y))); }

        public static List<T> LetReverse<T>(this List<T> list) { list.Reverse(); return list; }
        public static List<T> LetSort<T>(this List<T> list) where T : IComparable<T> { list.Sort(); return list; }
        public static List<T> LetSort<T>(this List<T> list, Comparison<T> compare) { list.Sort(compare); return list; }
        public static List<T> LetSort<T>(this List<T> list, Func<T, int> selector) { return list.LetSort((x, y) => selector(x).CompareTo(selector(y))); }
        public static List<T> LetSort<T>(this List<T> list, Func<T, double> selector) { return list.LetSort((x, y) => selector(x).CompareTo(selector(y))); }
        public static List<T> LetSort<T>(this List<T> list, Func<T, BigInteger> selector) { return list.LetSort((x, y) => selector(x).CompareTo(selector(y))); }

        public static int[] SortIndex<T>(this IList<T> list, Comparison<T> comparison)
        {
            int[] index = new int[list.Count];
            for (int i = index.Length; --i >= 0; )
                index[i] = i;
            Array.Sort(index, (x, y) => comparison(list[x], list[y]));
            return index;
        }
        public static int[] SortIndex<T>(this IList<T> list) where T : IComparable<T> { return SortIndex(list, (x, y) => x.CompareTo(y)); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, int> selector) { return SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, double> selector) { return SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))); }

        public static int[] IndexToRank(IList<int> index)
        {
            int[] Rank = new int[index.Count];
            for (int i = index.Count; --i >= 0; )
                Rank[index[i]] = i;
            return Rank;
        }

        public static T[] LetSortAsIndex<T>(this T[] list, IList<int> index) { return (T[])((IList<T>)list).LetSortAsIndex(index); }
        public static List<T> LetSortAsIndex<T>(this List<T> list, IList<int> index) { return (List<T>)((IList<T>)list).LetSortAsIndex(index); }
        public static IList<T> LetSortAsIndex<T>(this IList<T> list, IList<int> index)
        {
            for (int i = index.Count; --i >= 0; )
            {
                if (index[i] < 0) continue;
                T item = list[i];
                for (int j = i; ; )  // listの中の空いている添字番号
                {
                    int jj = index[j]; index[j] = ~index[j];
                    if (jj == i) { list[j] = item; break; }
                    list[j] = list[jj];
                    j = jj;
                }
            }
            for (int i = index.Count; --i >= 0; )
                index[i] = ~index[i];
            return list;
        }
        public static T[] OrderByIndex<T>(this IList<T> list, IList<int> index)
        {
            var result = new T[index.Count];
            for (int i = result.Length; --i >= 0; )
                result[i] = list[index[i]];
            return result;
        }

        public static T[] LetSortAsRank<T>(this T[] list, IList<int> rank) { return (T[])((IList<T>)list).LetSortAsRank(rank); }
        public static List<T> LetSortAsRank<T>(this List<T> list, IList<int> rank) { return (List<T>)((IList<T>)list).LetSortAsRank(rank); }
        public static IList<T> LetSortAsRank<T>(this IList<T> list, IList<int> rank)
        {
            for (int i = rank.Count; --i >= 0; )
            {
                if (rank[i] < 0) continue;
                T item = list[i];
                for (int j = i; ; )  // listの中の玉突きの玉の添字番号
                {
                    int jj = rank[j]; rank[j] = ~rank[j];
                    if (jj == i) { list[jj] = item; break; }
                    { T temp = list[jj]; list[jj] = item; item = temp; }
                    j = jj;
                }
            }
            for (int i = rank.Count; --i >= 0; )
                rank[i] = ~rank[i];
            return list;
        }
        public static T[] OrderByRank<T>(this IList<T> list, IList<int> rank)
        {
            var result = new T[list.Count];
            for (int i = result.Length; --i >= 0; )
                result[rank[i]] = list[i];
            return result;
        }

        public static int Compare<T>(IList<T> left, IList<T> right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return +1;
            if (left.Count != right.Count) return left.Count - right.Count;
            for (int i = 0; i < left.Count; i++)
            {
                int c = left[i].CompareTo(right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[] left, T[] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return +1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = left[i].CompareTo(right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[][] left, T[][] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return +1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = Compare(left[i], right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[][][] left, T[][][] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return +1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = Compare(left[i], right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[][][][] left, T[][][][] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return +1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = Compare(left[i], right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        #endregion

        #region TimeSpan
        public static TimeSpan Multiply(this TimeSpan left, int right) { return TimeSpan.FromTicks(left.Ticks * right); }
        public static TimeSpan Multiply(this TimeSpan left, double right) { return TimeSpan.FromTicks((long)(left.Ticks * right)); }
        public static TimeSpan Divide(this TimeSpan left, int right) { return TimeSpan.FromTicks(left.Ticks / right); }
        public static TimeSpan Divide(this TimeSpan left, double right) { return TimeSpan.FromTicks((long)(left.Ticks / right)); }
        public static TimeSpan Remainder(this TimeSpan left, int right) { return TimeSpan.FromTicks(left.Ticks % right); }
        public static TimeSpan Remainder(this TimeSpan left, double right) { return TimeSpan.FromTicks((long)(left.Ticks % right)); }
        #endregion

        #region string
        public static int TryParse(this string str, int defaultvalue)
        {
            int value;
            return int.TryParse(str, out value) ? value : defaultvalue;
        }
        public static double TryParse(this string str, double defaultvalue)
        {
            double value;
            return double.TryParse(str, out value) ? value : defaultvalue;
        }

        public static IEnumerable<string> SelectString<TSource>(this IEnumerable<TSource> source) { return source.Select(e => e.ToString()); }
        public static IEnumerable<string> SelectString<TSource>(this IEnumerable<TSource> source, string format)
        {
            var formatstr = "{0:" + format + "}";
            return source.Select(e => string.Format(formatstr, e));
        }

        public static string Join(this IEnumerable<string> source, string separator, string terminator = "")
        {
            var buffer = new StringBuilder();
            foreach (var element in source)
            {
                buffer.Append(element);
                buffer.Append(separator);
            }
            if (buffer.Length > 0) buffer.Remove(buffer.Length - separator.Length, separator.Length);
            buffer.Append(terminator);
            return buffer.ToString();
        }
        public static string JoinTab(this IEnumerable<string> source) { return source.Join("\t"); }
        public static string JoinSpace(this IEnumerable<string> source) { return source.Join(" "); }
        public static string JoinComma(this IEnumerable<string> source) { return source.Join(","); }
        public static string JoinCommaSpace(this IEnumerable<string> source) { return source.Join(", "); }
        public static string ToStringFormat<T>(this IEnumerable<T> source, string format = "", string separator = "\t", string terminator = "")
        {
            var formatstr = "{0:" + format + "}" + separator;
            var buffer = new StringBuilder();
            foreach (var element in source)
                buffer.AppendFormat(formatstr, element);
            if (buffer.Length > 0) buffer.Remove(buffer.Length - separator.Length, separator.Length);
            buffer.Append(terminator);
            return buffer.ToString();
        }
        public static string ToStringFormat<T>(this T[,] array, string format = "", string separator = "\t", string lineseparator = "\n", string terminator = "")
        {
            var formatstr = "{0:" + format + "}" + separator;
            var buffer = new StringBuilder();
            for (int i0 = 0; i0 < array.GetLength(0); i0++)
            {
                for (int i1 = 0; i1 < array.GetLength(1); i1++)
                    buffer.AppendFormat(formatstr, array[i0, i1]);
                buffer.Remove(buffer.Length - separator.Length, separator.Length);
                buffer.Append(lineseparator);
            }
            buffer.Append(terminator);
            return buffer.ToString();
        }
        public static string ToStringFormat<T>(this T[, ,] array, string format = "", string separator = "\t", string lineseparator = "\n", string terminator = "")
        {
            var formatstr = "{0:" + format + "}" + separator;
            var buffer = new StringBuilder();
            for (int i0 = 0; i0 < array.GetLength(0); i0++)
            {
                for (int i1 = 0; i1 < array.GetLength(1); i1++)
                {
                    for (int i2 = 0; i2 < array.GetLength(1); i2++)
                        buffer.AppendFormat(formatstr, array[i0, i1, i2]);
                    buffer.Remove(buffer.Length - separator.Length, separator.Length);
                    buffer.Append(lineseparator);
                }
            }
            buffer.Append(terminator);
            return buffer.ToString();
        }

        public static string[] Split(this string str, char separator) { return str.Split(new char[] { separator }); }
        public static string[] Split(this string str, string separator) { return str.Split(new string[] { separator }, StringSplitOptions.None); }
        public static string[] SplitLine(this string str) { return str.Split(Ct.NewLineCodes, StringSplitOptions.None); }
        public static string[] SplitTab(this string line) { return line.Split('\t'); }
        public static string[] SplitCsv(this string line)
        {
            if (line == null) return null;
            var converted = new List<char>();
            bool quotation = false;
            bool escape = false;
            for (int i = 0; i < line.Length; i++)
            {
                char c = i < line.Length ? line[i] : ',';
                if (escape)
                {
                    if (c != '\"') { quotation = false; escape = false; }
                }
                if (!quotation)
                {
                    if (c == '\"') { quotation = true; continue; }
                    if (c == ',') c = '\t';
                }
                else
                {
                    if (c == '\"')
                    {
                        if (!escape) { escape = true; continue; }  // 1個目
                        else { escape = false; }                   // 2個目
                    }
                }
                converted.Add(c);
            }
            return new string(converted.ToArray()).SplitTab();
        }

        public static string Replace(this string str, string[,] replacelist)
        {
            for (int i = 0; i < replacelist.GetLength(0); i++)
                str = str.Replace(replacelist[i, 0], replacelist[i, 1]);
            return str;
        }

        public static void Write(this StringBuilder sb, string value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, bool value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, char value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, char[] buffer) { sb.Append(buffer, 0, buffer.Length); }
        public static void Write(this StringBuilder sb, char[] buffer, int index, int count) { sb.Append(buffer, index, count); }
        public static void Write(this StringBuilder sb, decimal value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, double value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, float value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, int value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, long value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, object value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, string format, object arg0) { sb.Append(string.Format(format, arg0)); }
        public static void Write(this StringBuilder sb, string format, object arg0, object arg1) { sb.Append(string.Format(format, arg0, arg1)); }
        public static void Write(this StringBuilder sb, string format, object arg0, object arg1, object arg2) { sb.Append(string.Format(format, arg0, arg1, arg2)); }
        public static void Write(this StringBuilder sb, string format, params object[] arg) { sb.Append(string.Format(format, arg)); }
        public static void Write(this StringBuilder sb, uint value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, ulong value) { sb.Append(value); }

        public static void WriteLine(this StringBuilder sb) { sb.AppendLine(); }
        public static void WriteLine(this StringBuilder sb, string value) { sb.AppendLine(value); }
        public static void WriteLine(this StringBuilder sb, bool value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, char value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, char[] buffer) { sb.Append(buffer); sb.AppendLine(); }
        public static void WriteLine(this StringBuilder sb, char[] buffer, int index, int count) { sb.Append(buffer, index, count); sb.AppendLine(); }
        public static void WriteLine(this StringBuilder sb, decimal value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, double value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, float value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, int value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, long value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, object value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, string format, object arg0) { sb.AppendLine(string.Format(format, arg0)); }
        public static void WriteLine(this StringBuilder sb, string format, object arg0, object arg1) { sb.AppendLine(string.Format(format, arg0, arg1)); }
        public static void WriteLine(this StringBuilder sb, string format, object arg0, object arg1, object arg2) { sb.AppendLine(string.Format(format, arg0, arg1, arg2)); }
        public static void WriteLine(this StringBuilder sb, string format, params object[] arg) { sb.AppendLine(string.Format(format, arg)); }
        public static void WriteLine(this StringBuilder sb, uint value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, ulong value) { sb.AppendLine(value.ToString()); }
        #endregion

        #region IO
        public static string[] GetFiles(string path) { return GetFiles(path, SearchOption.TopDirectoryOnly); }
        public static string[] GetFiles(string path, SearchOption option)
        {
            string dirname = Path.GetDirectoryName(path);
            string filename = Path.GetFileName(path);
            if (filename.Length == 0) filename = "*";
            string[] paths = Directory.GetFiles(dirname, filename, option);
            return paths;
        }

        public static void Write(this Stream stream, IEnumerable<byte> source) { stream.Write(source.ToArray()); }
        public static void Write(this Stream stream, byte[] array) { stream.Write(array, 0, array.Length); }
        public static int Read(this FileStream filestream, byte[] array) { return filestream.Read(array, 0, array.Length); }
        public static bool EndOfStream(this Stream stream) { return stream.Position >= stream.Length; }
        public static bool EndOfStream(this BinaryReader reader) { return reader.BaseStream.EndOfStream(); }

        #endregion
    }

    // Unsafe array functions
    public unsafe static partial class Us
    {
        #region double
        public static double Sum(double* p, int n, Func<double, double> f) { double a = 0; for (int i = n; --i >= 0; ) a += f(p[i]); return a; }
        public static double Sum(double* p, double* q, int n, Func<double, double, double> f) { double a = 0; for (int i = n; --i >= 0; ) a += f(p[i], q[i]); return a; }
        public static double Sum(double* p, int n) { double a = 0; for (int i = n; --i >= 0; ) a += p[i]; return a; }
        public static double SumAbs(double* p, int n) { double a = 0; for (int i = n; --i >= 0; ) a += Math.Abs(p[i]); return a; }
        public static double SumSq(double* p, int n) { double a = 0; for (int i = n; --i >= 0; ) a += p[i] * p[i]; return a; }
        public static double SumPow(double* p, int n, double nu) { double a = 0; for (int i = n; --i >= 0; ) a += Math.Pow(Math.Abs(p[i]), nu); return a; }
        public static double SumMul(double* p, double* q, int n) { double a = 0; for (int i = n; --i >= 0; ) a += p[i] * q[i]; return a; }
        public static double SumAbsSub(double* p, double* q, int n) { double a = 0; for (int i = n; --i >= 0; ) a += Math.Abs(p[i] - q[i]); return a; }
        public static double SumSqSub(double* p, double* q, int n) { double a = 0; for (int i = n; --i >= 0; ) a += Mt.Sq(p[i] - q[i]); return a; }
        public static double SumPowSub(double* p, double* q, int n, double nu) { double a = 0; for (int i = n; --i >= 0; ) a += Math.Pow(Math.Abs(p[i] - q[i]), nu); return a; }
        public static double Min(double* p, int n, Func<double, double> f) { double a = 0; for (int i = 0; i < n; i++) { var b = f(p[i]); if (i == 0 || a > b) a = b; } return a; }
        public static double Max(double* p, int n, Func<double, double> f) { double a = 0; for (int i = 0; i < n; i++) { var b = f(p[i]); if (i == 0 || a < b) a = b; } return a; }
        public static double Min(double* p, double* q, int n, Func<double, double, double> f) { double a = 0; for (int i = 0; i < n; i++) { var b = f(p[i], q[i]); if (i == 0 || a > b) a = b; } return a; }
        public static double Max(double* p, double* q, int n, Func<double, double, double> f) { double a = 0; for (int i = 0; i < n; i++) { var b = f(p[i], q[i]); if (i == 0 || a < b) a = b; } return a; }
        public static double Min(double* p, int n) { double a = 0; for (int i = 0; i < n; i++) { var b = p[i]; if (i == 0 || a > b) a = b; } return a; }
        public static double Max(double* p, int n) { double a = 0; for (int i = 0; i < n; i++) { var b = p[i]; if (i == 0 || a < b) a = b; } return a; }
        public static double MinAbs(double* p, int n) { double a = 0; for (int i = 0; i < n; i++) { var b = Math.Abs(p[i]); if (i == 0 || a > b) a = b; } return a; }
        public static double MaxAbs(double* p, int n) { double a = 0; for (int i = 0; i < n; i++) { var b = Math.Abs(p[i]); if (i == 0 || a < b) a = b; } return a; }
        public static double MinSq(double* p, int n) { double a = 0; for (int i = 0; i < n; i++) { var b = Mt.Sq(p[i]); if (i == 0 || a > b) a = b; } return a; }
        public static double MaxSq(double* p, int n) { double a = 0; for (int i = 0; i < n; i++) { var b = Mt.Sq(p[i]); if (i == 0 || a < b) a = b; } return a; }
        public static double MinPow(double* p, int n, double nu) { double a = 0; for (int i = 0; i < n; i++) { var b = Math.Pow(Math.Abs(p[i]), nu); if (i == 0 || a > b) a = b; } return a; }
        public static double MaxPow(double* p, int n, double nu) { double a = 0; for (int i = 0; i < n; i++) { var b = Math.Pow(Math.Abs(p[i]), nu); if (i == 0 || a < b) a = b; } return a; }

        public static void LetStep(double* r, double* p, int step, int n) { for (int i = 0, j = 0; i < n; i++, j += step) r[i] = p[j]; }
        public static void LetStep(double* r, int step, double* p, int n) { for (int i = 0, j = 0; i < n; i++, j += step) r[j] = p[i]; }
        public static void LetStartStep(double* r, double* p, int start, int step, int n) { for (int i = 0, j = start; i < n; i++, j += step) r[i] = p[j]; }
        public static void LetStartStep(double* r, int start, int step, double* p, int n) { for (int i = 0, j = start; i < n; i++, j += step) r[j] = p[i]; }
        public static void Let(double* r, double* p, int n) { for (int i = n; --i >= 0; ) r[i] = p[i]; }
        public static void Let(double* r, double p, int n) { for (int i = n; --i >= 0; ) r[i] = p; }
        public static void LetRev(double* r, int n) { for (int i = n / 2; --i >= 0; ) { var a = r[n - 1 - i]; r[n - 1 - i] = r[i]; r[i] = a; } }
        public static void Neg(double* r, double* p, int n) { for (int i = n; --i >= 0; ) r[i] = -p[i]; }
        public static void Rev(double* r, double* p, int n) { for (int i = n; --i >= 0; ) r[i] = p[n - 1 - i]; }
        public static void Add(double* r, double* p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] + q[i]; }
        public static void Sub(double* r, double* p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] - q[i]; }
        public static void Mul(double* r, double* p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] * q[i]; }
        public static void Div(double* r, double* p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] / q[i]; }
        public static void Mod(double* r, double* p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] % q[i]; }
        public static void Add(double* r, double* p, double q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] + q; }
        public static void Sub(double* r, double* p, double q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] - q; }
        public static void Mul(double* r, double* p, double q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] * q; }
        public static void Div(double* r, double* p, double q, int n) { Mul(r, p, 1 / q, n); }
        public static void Mod(double* r, double* p, double q, int n) { for (int i = n; --i >= 0; ) r[i] = p[i] % q; }
        public static void Sub(double* r, double p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p - q[i]; }
        public static void Div(double* r, double p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p / q[i]; }
        public static void Mod(double* r, double p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = p % q[i]; }

        public static void Operate(double* r, double* p, int n, Func<double, double> f) { for (int i = n; --i >= 0; ) r[i] = f(p[i]); }
        public static void Operate(double* r, double* p, double* q, int n, Func<double, double, double> f) { for (int i = n; --i >= 0; ) r[i] = f(p[i], q[i]); }

        public static void LetAddMul(double* r, double* p, double q, int n)
        {
            int m = n & ~3;
            for (int i = n; --i >= m; )
                r[i] += p[i] * q;
            for (int i = m; (i -= 4) >= 0; )
            {
                *(r + i + 3) += *(p + i + 3) * q;
                *(r + i + 2) += *(p + i + 2) * q;
                *(r + i + 1) += *(p + i + 1) * q;
                *(r + i + 0) += *(p + i + 0) * q;
            }
        }
        public static void LetMulAdd(double* r, double p, double* q, int n) { for (int i = n; --i >= 0; ) r[i] = r[i] * p + q[i]; }
        public static void LetMulAdd(double* r, double p, double q, int n) { for (int i = n; --i >= 0; ) r[i] = r[i] * p + q; }
        #endregion

        #region complex
        public static double SqAbs(complex* p) { return p->Re * p->Re + p->Im * p->Im; }
        public static void Swap(complex* p, complex* q)
        {
            double tr = p->Re;
            double ti = p->Im;
            p->Re = q->Re;
            p->Im = q->Im;
            q->Re = tr;
            q->Im = ti;
        }
        public static void Let(complex* p, complex* q)
        {
            p->Re = q->Re;
            p->Im = q->Im;
        }
        public static void Div(complex* r, double p, complex* q)
        {
            double tl = p / (q->Re * q->Re + q->Im * q->Im);
            r->Re = +q->Re * tl;
            r->Im = -q->Im * tl;
        }
        public static void LetAddMul(complex* r, complex* p, complex* q)
        {
            double tr = p->Re * q->Re - p->Im * q->Im;
            double ti = p->Re * q->Im + p->Im * q->Re;
            r->Re += tr;
            r->Im += ti;
        }
        public static void LetAddCul(complex* r, complex* p, complex* q)
        {
            double tr = p->Re * q->Re + p->Im * q->Im;
            double ti = p->Re * q->Im - p->Im * q->Re;
            r->Re += tr;
            r->Im += ti;
        }
        public static void Mul(complex* r, complex* p, complex* q)
        {
            double tr = p->Re * q->Re - p->Im * q->Im;
            double ti = p->Re * q->Im + p->Im * q->Re;
            r->Re = tr;
            r->Im = ti;
        }
        public static void Div(complex* r, complex* p, complex* q)
        {
            double tl = q->Re * q->Re + q->Im * q->Im;
            double tr = (p->Im * q->Im + p->Re * q->Re) / tl;
            double ti = (p->Im * q->Re - p->Re * q->Im) / tl;
            r->Re = tr;
            r->Im = ti;
        }
        public static void Cul(complex* r, complex* p, complex* q)
        {
            double tr = p->Re * q->Re + p->Im * q->Im;
            double ti = p->Re * q->Im - p->Im * q->Re;
            r->Re = tr;
            r->Im = ti;
        }
        public static double SqAbsSub(complex* p, complex* q) { return Mt.Sq(p->Re - q->Re) + Mt.Sq(p->Im - q->Im); }

        public static double SumAbs(complex* p, int n) { double a = 0; for (int i = n; --i >= 0; ) a += Math.Sqrt(SqAbs(&p[i])); return a; }
        public static double SumSqAbs(complex* p, int n) { double a = 0; for (int i = n; --i >= 0; ) a += SqAbs(&p[i]); return a; }
        public static double MaxSqAbs(complex* p, int n) { double a = 0; for (int i = n; --i >= 0; ) Mt.LetMax(ref a, SqAbs(&p[i])); return a; }
        public static double SumPowAbs(complex* p, int n, double nu) { nu *= 0.5; double a = 0; for (int i = n; --i >= 0; ) a += Math.Pow(SqAbs(&p[i]), nu); return a; }
        public static complex Sum(complex* p, int n) { complex a = default(complex); for (int i = n; --i >= 0; ) { a.Im += p[i].Im; a.Re += p[i].Re; } return a; }
        public static complex SumMul(complex* p, complex* q, int n) { var a = default(complex); for (int i = n; --i >= 0; ) LetAddMul(&a, &p[i], &q[i]); return a; }
        public static complex SumCul(complex* p, complex* q, int n) { var a = default(complex); for (int i = n; --i >= 0; ) LetAddCul(&a, &p[i], &q[i]); return a; }
        public static complex SumMul(complex* p, double* q, int n) { var a = default(complex); for (int i = n; --i >= 0; ) { a.Im += p[i].Im * q[i]; a.Re += p[i].Re * q[i]; } return a; }
        public static complex SumCul(complex* p, double* q, int n) { var a = SumMul(p, q, n); a.Im *= -1; return a; }
        public static double SumAbsSub(complex* p, complex* q, int n) { double a = 0; for (int i = n; --i >= 0; ) a += Math.Sqrt(SqAbsSub(&p[i], &q[i])); return a; }
        public static double SumSqSub(complex* p, complex* q, int n) { double a = 0; for (int i = n; --i >= 0; ) a += SqAbsSub(&p[i], &q[i]); return a; }
        public static double SumPowSub(complex* p, complex* q, int n, double nu) { nu *= 0.5; double a = 0; for (int i = n; --i >= 0; ) a += Math.Pow(SqAbsSub(&p[i], &q[i]), nu); return a; }

        public static void Tocomplex(complex* r, double* p, double* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = q[i]; r[i].Re = p[i]; } }
        public static void Real(double* r, complex* p, int n) { for (int i = n; --i >= 0; ) r[i] = p[i].Re; }
        public static void Imag(double* r, complex* p, int n) { for (int i = n; --i >= 0; ) r[i] = p[i].Im; }
        public static void LetStartStep(complex* r, complex* p, int start, int step, int n) { for (int i = 0, j = start; i < n; i++, j += step) { r[i].Re = p[j].Re; r[i].Im = p[j].Im; } }
        public static void LetStartStep(complex* r, int start, int step, complex* p, int n) { for (int i = 0, j = start; i < n; i++, j += step) { r[j].Re = p[i].Re; r[j].Im = p[i].Im; } }
        public static void Let(complex* p, complex* q, int n) { Let((double*)p, (double*)q, n * 2); }
        public static void Let(complex* p, double* q, int n) { for (int i = n; --i >= 0; ) { p[i].Im = 0; p[i].Re = q[i]; } }
        public static void Let(complex* r, complex p, int n) { for (int i = n; --i >= 0; ) r[i] = p; }
        public static void LetRev(complex* r, int n) { for (int i = n / 2; --i >= 0; ) { var a = r[n - 1 - i]; r[n - 1 - i] = r[i]; r[i] = a; } }
        public static void Neg(complex* r, complex* p, int n) { Neg((double*)r, (double*)p, n * 2); }
        public static void Rev(complex* r, complex* p, int n) { for (int i = n; --i >= 0; ) r[i] = p[n - 1 - i]; }
        public static void Conj(complex* r, complex* p, int n) { for (int i = n; --i >= 0; ) { r[i].Im = -p[i].Im; r[i].Re = p[i].Re; } }
        public static void Add(complex* r, complex* p, complex* q, int n) { Add((double*)r, (double*)p, (double*)q, n * 2); }
        public static void Sub(complex* r, complex* p, complex* q, int n) { Sub((double*)r, (double*)p, (double*)q, n * 2); }
        public static void Mul(complex* r, complex* p, complex* q, int n) { for (int i = n; --i >= 0; ) Mul(&r[i], &p[i], &q[i]); }
        public static void Div(complex* r, complex* p, complex* q, int n) { for (int i = n; --i >= 0; ) Div(&r[i], &p[i], &q[i]); }
        public static void Cul(complex* r, complex* p, complex* q, int n) { for (int i = n; --i >= 0; ) Cul(&r[i], &p[i], &q[i]); }
        public static void Add(complex* r, complex* p, complex q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im + q.Im; r[i].Re = p[i].Re + q.Re; } }
        public static void Sub(complex* r, complex* p, complex q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im - q.Im; r[i].Re = p[i].Re - q.Re; } }
        public static void Mul(complex* r, complex* p, complex q, int n) { for (int i = n; --i >= 0; ) Mul(&r[i], &p[i], &q); }
        public static void Div(complex* r, complex* p, complex q, int n) { Div(&q, 1, &q); Mul(r, p, q, n); }
        public static void Cul(complex* r, complex* p, complex q, int n) { for (int i = n; --i >= 0; ) Cul(&r[i], &p[i], &q); }
        public static void Sub(complex* r, complex p, complex* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p.Im - q[i].Im; r[i].Re = p.Re - q[i].Re; } }
        public static void Div(complex* r, complex p, complex* q, int n) { for (int i = n; --i >= 0; ) Div(&r[i], &p, &q[i]); }

        public static void Add(complex* r, complex* p, double* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im; r[i].Re = p[i].Re + q[i]; } }
        public static void Sub(complex* r, complex* p, double* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im; r[i].Re = p[i].Re - q[i]; } }
        public static void Mul(complex* r, complex* p, double* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im * q[i]; r[i].Re = p[i].Re * q[i]; } }
        public static void Div(complex* r, complex* p, double* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im / q[i]; r[i].Re = p[i].Re / q[i]; } }
        public static void Cul(complex* r, complex* p, double* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = -p[i].Im * q[i]; r[i].Re = p[i].Re * q[i]; } }
        public static void Sub(complex* r, double* p, complex* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = -q[i].Im; r[i].Re = p[i] - q[i].Re; } }
        public static void Div(complex* r, double* p, complex* q, int n) { for (int i = n; --i >= 0; ) Div(&r[i], p[i], &q[i]); }
        public static void Add(complex* r, complex* p, double q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im; r[i].Re = p[i].Re + q; } }
        public static void Sub(complex* r, complex* p, double q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = p[i].Im; r[i].Re = p[i].Re - q; } }
        public static void Mul(complex* r, complex* p, double q, int n) { Mul((double*)r, (double*)p, q, n * 2); }
        public static void Div(complex* r, complex* p, double q, int n) { Div((double*)r, (double*)p, q, n * 2); }
        public static void Cul(complex* r, complex* p, double q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = -p[i].Im * q; r[i].Re = p[i].Re * q; } }
        public static void Sub(complex* r, double p, complex* q, int n) { for (int i = n; --i >= 0; ) { r[i].Im = -q[i].Im; r[i].Re = p - q[i].Re; } }
        public static void Div(complex* r, double p, complex* q, int n) { for (int i = n; --i >= 0; ) Div(&r[i], p, &q[i]); }

        public static void LetAddMul(complex* r, complex* p, complex q, int n) { for (int i = n; --i >= 0; ) LetAddMul(&r[i], &p[i], &q); }
        public static void LetAddMul(complex* r, complex* p, double q, int n) { LetAddMul((double*)r, (double*)p, q, n * 2); }
        public static void LetMulAdd(complex* r, double p, complex* q, int n) { LetMulAdd((double*)r, p, (double*)q, n * 2); }
        public static void SqAbs(double* r, complex* p, int n) { for (int i = n; --i >= 0; ) r[i] = SqAbs(&p[i]); }
        #endregion
    }

    // Mathematical functions
    public static partial class Mt
    {
        #region Miscellaneous functions
        public const double DoubleEpsilon = 2.2250738585072013830902327173324e-308;  // 2^-1022
        public const double DoubleEps = 2.2204460492503130808472633361816e-16;  // 2^-52
        public const double DoubleFpMin = 1.0020841800044863889980540256751e-292;  // 2^-1022 / 2^-52 = 2^-970
        public const double PI2 = Math.PI * 2;
        public const double Ln2 = 0.69314718055994530941723212145818; // Math.Log(2);
        public const double Sqrt2 = 1.4142135623730950488016887242097; // Math.Sqrt(2);
        public const double Sqrt2inv = 0.70710678118654752440084436210485;  // 1 / Math.Sqrt(2);
        public const double Sqrt2PIinv = 0.39894228040143267793994605993438;  // 1 / Math.Sqrt(2 * Math.PI);

        public static bool IsNaN(this double value) { return double.IsNaN(value); }
        public static bool IsInf(this double value) { return double.IsInfinity(value); }
        public static bool IsPosInf(this double value) { return double.IsPositiveInfinity(value); }
        public static bool IsNegInf(this double value) { return double.IsNegativeInfinity(value); }
        public static bool IsNanOrInf(this double value) { return double.IsInfinity(value) || double.IsNaN(value); }
        public static bool ContainsNanOrInf(this IEnumerable<double> source) { return source.Any(v => v.IsNanOrInf()); }

        public static bool IsTooSmall(this double x, double y) { return x + y == y; }

        public static void Swap<T>(ref T val1, ref T val2) { T z = val1; val1 = val2; val2 = z; }

        public static T MinMax<T>(T value, T valmin, T valmax) where T : IComparable<T> { return value.CompareTo(valmin) < 0 ? valmin : (value.CompareTo(valmax) > 0 ? valmax : value); }
        public static void LetMin<T>(ref T value, T valmin) where T : IComparable<T> { if (value.CompareTo(valmin) > 0) value = valmin; }
        public static void LetMax<T>(ref T value, T valmax) where T : IComparable<T> { if (value.CompareTo(valmax) < 0) value = valmax; }
        public static void LetMinMax<T>(ref T value, T valmin, T valmax) where T : IComparable<T> { if (value.CompareTo(valmin) < 0) value = valmin; else if (value.CompareTo(valmax) > 0) value = valmax; }
        public static bool IfLetMin<T>(ref T value, T valmin) where T : IComparable<T> { if (value.CompareTo(valmin) > 0) { value = valmin; return true; } return false; }
        public static bool IfLetMax<T>(ref T value, T valmax) where T : IComparable<T> { if (value.CompareTo(valmax) < 0) { value = valmax; return true; } return false; }
        public static int MinMaxC(int value, int valmin, int valmax) { return value < valmin ? valmin : value >= valmax ? valmax - 1 : value; }
        public static void LetMinMaxC(ref int value, int valmin, int valmax) { if (value < valmin) value = valmin; else if (value >= valmax) value = valmax - 1; }

        public static Int2 DivRem(int left, int right) { int rem, div = Math.DivRem(left, right, out rem); return new Int2(rem, div); }
        public static int AlignUp(int value, int size) { return (value + (size - 1)) - (value + (size - 1)) % size; }
        public static int AlignDown(int value, int size) { return value - value % size; }
        public static int AlignUpX(int value, int size) { return (value + (size - 1)) & (~(size - 1)); }
        public static int AlignDownX(int value, int size) { return value & (~(size - 1)); }
        public static int RoundOff(double value) { return (int)Math.Floor(value + 0.5); }
        public static double RoundOff(double value, int digit) { double k = Math.Pow(10, digit); return Math.Floor(value * k + 0.5) / k; }
        public static Int2 RoundOff(Double2 value) { return new Int2(Mt.RoundOff(value.X), Mt.RoundOff(value.Y)); }
        public static Int3 RoundOff(Double3 value) { return new Int3(Mt.RoundOff(value.X), Mt.RoundOff(value.Y), Mt.RoundOff(value.Z)); }

        public static int Div_(int left, int right) { return (left - ((left % right) + right) % right) / right; }  //right>0:下方向商   right<0:上方向商
        public static int Mod_(int left, int right) { return ((left % right) + right) % right; }        //right>0:下方向剰余 right<0:上方向剰余        
        public static double Div_(double left, double right) { return (left - ((left % right) + right) % right) / right; }  //right>0:下方向商   right<0:上方向商
        public static double Mod_(double left, double right) { return ((left % right) + right) % right; }        //right>0:下方向剰余 right<0:上方向剰余        
        //  5/ 3= 1	Div_( 5, 3)= 1
        // -5/ 3=-1	Div_(-5, 3)=-2
        //  5/-3=-1	Div_( 5,-3)=-2
        // -5/-3= 1	Div_(-5,-3)= 1
        //  5% 3= 2	Mod_( 5, 3)= 2
        // -5% 3=-2	Mod_(-5, 3)= 1
        //  5%-3= 2	Mod_( 5,-3)=-1
        // -5%-3=-2	Mod_(-5,-3)=-2
        public static double ModPmHalf(double value)
        {
            return Mod_((value - 0.5), -1) + 0.5;
        }
        public static double ModPmPi(double value)
        {
            return Mod_((value - Math.PI), -PI2) + Math.PI;
        }

        public static int Sq(this int value) { return value * value; }
        public static double Sq(this double value) { return value * value; }
        public static int Cube(this int value) { return value * value * value; }
        public static double Cube(this double value) { return value * value * value; }
        public static int LetAbs(ref int value) { if (value < 0) value *= -1; return value; }
        public static double LetAbs(ref double value) { if (value < 0) value *= -1; return value; }
        public static void LetOrder(ref int val1, ref int val2) { if (val2 < val1) Swap(ref val1, ref val2); }
        public static void LetOrder(ref double val1, ref double val2) { if (val2 < val1) Swap(ref val1, ref val2); }

        public static int DivCeil(int left, int right) { return (left + right - 1) / right; }
        public static long DivCeil(long left, int right) { return (left + right - 1) / right; }
        public static int DivRound(int left, int right) { return (left + (left >= 0 ? right / 2 : -right / 2)) / right; }
        public static long DivRound(long left, int right) { return (left + (left >= 0 ? right / 2 : -right / 2)) / right; }

        public static int Min(int val1, int val2, int val3) { return Math.Min(Math.Min(val1, val2), val3); }
        public static int Max(int val1, int val2, int val3) { return Math.Max(Math.Max(val1, val2), val3); }
        public static int Min(int val1, int val2, int val3, int val4) { return Math.Min(Math.Min(val1, val2), Math.Min(val3, val4)); }
        public static int Max(int val1, int val2, int val3, int val4) { return Math.Max(Math.Max(val1, val2), Math.Max(val3, val4)); }
        public static double Min(double val1, double val2, double val3) { return Math.Min(Math.Min(val1, val2), val3); }
        public static double Max(double val1, double val2, double val3) { return Math.Max(Math.Max(val1, val2), val3); }
        public static double Min(double val1, double val2, double val3, double val4) { return Math.Min(Math.Min(val1, val2), Math.Min(val3, val4)); }
        public static double Max(double val1, double val2, double val3, double val4) { return Math.Max(Math.Max(val1, val2), Math.Max(val3, val4)); }

        public static double SlightlyInferior(double value) { return value - Math.Max(1e-296, Math.Abs(value) * 1e-12); }
        public static double SlightlySuperior(double value) { return value + Math.Max(1e-296, Math.Abs(value) * 1e-12); }
        public static int Log2Floor(int value)
        {
            if (value <= 0) ThrowException.ArgumentOutOfRange("value");
            return Log2Floor((uint)value);
        }
        public static int Log2Floor(uint value)
        {
            if (value == 0) ThrowException.ArgumentOutOfRange("value");
            int i = 16;
            for (int j = 8; j > 0; j >>= 1)
                if (value < (1u << i)) i -= j; else i += j;
            return value < (1u << i) ? i - 1 : i;
        }
        public static int Log2Ceiling(int value)
        {
            if (value <= 0) ThrowException.ArgumentOutOfRange("value");
            return Log2Ceiling((uint)value);
        }
        public static int Log2Ceiling(uint value)
        {
            if (value == 0) ThrowException.ArgumentOutOfRange("value");
            if (value == 1) return 0;
            int i = 16;
            for (int j = 8; j > 0; j >>= 1)
                if (value <= (1u << i)) i -= j; else i += j;
            return value > (1u << i) ? i + 1 : i;
        }

        public static double WLog(double weight, double value)
        {
            if (value == 0 && weight == 0) return 0;
            return weight * Math.Log(value);
        }
        public static double Logistic(double value)
        {
            return 1 / (1 + Math.Exp(-value));
        }
        public static double Atanh(double value)
        {
            return (Math.Abs(value) >= 1 ? (value > 0 ? double.PositiveInfinity : double.NegativeInfinity) : 0.5 * Math.Log((1 + value) / (1 - value)));
        }
        public static double Atanh_(double value)
        {
            return (Math.Abs(value) >= 1 ? (value > 0 ? 18.7149738751185 : -18.7149738751185) : 0.5 * Math.Log((1 + value) / (1 - value)));
        }
        public static double Tanh(double value)
        {
            if (value == double.PositiveInfinity) return +1;
            if (value == double.NegativeInfinity) return -1;
            double y = Math.Exp(value * 2);
            return (y - 1) / (y + 1);
        }
        public static double Tanh_(double value)
        {
            if (Math.Abs(value) >= 18.7149738751185) return value > 0 ? +1 - 1e-16 : -1 + 1e-16;
            double y = Math.Exp(value + value);
            return (y - 1) / (y + 1);
        }
        public static double Acos_(double value)
        {
            return value <= -1 ? Math.PI : value >= 1 ? 0.0 : Math.Acos(value);
        }
        public static double Asin_(double value)
        {
            return value <= -1 ? -0.5 * Math.PI : value >= 1 ? 0.5 * Math.PI : Math.Asin(value);
        }
        public static double Atan1(double y, double x)
        {
            double theta = Math.Atan2(y, x);
            if (theta <= -0.5 * Math.PI) return theta + Math.PI;
            if (theta > 0.5 * Math.PI) return theta - Math.PI;
            return theta;
        }

        //Sqrt(x*x + y*y) をoverflow, underflowに気をつけて計算
        public static double Norm2(double x, double y)
        {
            return double.IsInfinity(x) || double.IsInfinity(y) ? double.PositiveInfinity :
                Math.Abs(x) > Math.Abs(y) ?
                Math.Sqrt(1 + Sq(y / x)) * Math.Abs(x) :
                y == 0 ? Math.Abs(x) :
                Math.Sqrt(1 + Sq(x / y)) * Math.Abs(y);
        }

        public static complex Phase(double phase) { return new complex(Math.Cos(phase), Math.Sin(phase)); }
        // denominator > 0
        public static complex Phase(int numerator, int denominator)
        {
            numerator %= denominator;
            if (numerator >= 0)
            {
                if (numerator * 2 > denominator) numerator -= denominator;
            }
            else
            {
                if (numerator * -2 > denominator) numerator += denominator;
            }
            return Phase(PI2 * numerator / denominator);
        }
        public static double SqAbs(this complex value) { return value.Re * value.Re + value.Im * value.Im; }
        public static double Abs(this complex value) { return Math.Sqrt(value.Re * value.Re + value.Im * value.Im); }
        public static double Arg(this complex value) { return Math.Atan2(value.Im, value.Re); }
        public static complex Conj(this complex value) { return new complex(value.Re, -value.Im); }
        public static complex ConjMul(complex value0, complex value1)
        {
            return new complex(
                value0.Re * value1.Re + value0.Im * value1.Im,
                value0.Re * value1.Im - value0.Im * value1.Re);
        }
        public static complex ConjMul(complex value0, double value1) { return new complex(value0.Re * value1, -value0.Im * value1); }

        #endregion

        #region Basic calculations
        #region Int32
        public static double Average(int count, Func<int, int> function) { return (double)Sum(count, function) / count; }
        public static int Sum(int count, Func<int, int> function)
        {
            var a = default(int);
            for (int i = 0; i < count; i++) a += function(i);
            return a;
        }
        public static double GeometricalAverage(int count, Func<int, int> function) { return Math.Pow((double)Product(count, function), 1.0 / count); }
        public static int Product(int count, Func<int, int> function)
        {
            int a = 1;
            for (int i = 0; i < count; i++) checked { a *= function(i); }
            return a;
        }
        #endregion

        #region BigInteger
        public static double Average(int count, Func<int, BigInteger> function) { return (double)Sum(count, function) / count; }
        public static BigInteger Sum(int count, Func<int, BigInteger> function)
        {
            var a = default(BigInteger);
            for (int i = 0; i < count; i++) a += function(i);
            return a;
        }
        public static double GeometricalAverage(int count, Func<int, BigInteger> function) { return Math.Pow((double)Product(count, function), 1.0 / count); }
        public static BigInteger Product(int count, Func<int, BigInteger> function)
        {
            BigInteger a = 1;
            for (int i = 0; i < count; i++) a *= function(i);
            return a;
        }
        #endregion

        #region Double
        public static double Average(int count, Func<int, double> function) { return Sum(count, function) / count; }
        public static double Sum(int count, Func<int, double> function)
        {
            var a = default(double);
            for (int i = 0; i < count; i++) a += function(i);
            return a;
        }
        public static double GeometricalAverage(int count, Func<int, double> function) { return Math.Pow(Product(count, function), 1.0 / count); }
        public static double Product(int count, Func<int, double> function)
        {
            double a = 1;
            for (int i = 0; i < count; i++) a *= function(i);
            return a;
        }
        #endregion

        #region Double[]
        public static double[] Average(int count, Func<int, double[]> function) { return Sum(count, function).LetDiv(count); }
        public static double[] Sum(int count, Func<int, double[]> function)
        {
            var a = (count == 0) ? default(double[]) : function(0).Pos();
            for (int i = 1; i < count; i++) a.LetAdd(function(i));
            return a;
        }
        public static double[] Sum(int count, double[] init, Func<int, double[]> function)
        {
            var a = init;
            for (int i = 0; i < count; i++) a.LetAdd(function(i));
            return a;
        }
        #endregion

        #region Double[,]
        public static double[,] Average(int count, Func<int, double[,]> function) { return Sum(count, function).LetDiv(count); }
        public static double[,] Sum(int count, Func<int, double[,]> function)
        {
            var a = (count == 0) ? default(double[,]) : function(0).Pos();
            for (int i = 1; i < count; i++) a.LetAdd(function(i));
            return a;
        }
        #endregion

        #region complex
        public static complex Average(int count, Func<int, complex> function) { return Sum(count, function) / count; }
        public static complex Sum(int count, Func<int, complex> function)
        {
            var a = default(complex);
            for (int i = 0; i < count; i++) a += function(i);
            return a;
        }
        public static complex GeometricalAverage(int count, Func<int, complex> function) { return complex.Pow(Product(count, function), 1.0 / count); }
        public static complex Product(int count, Func<int, complex> function)
        {
            complex a = complex.One;
            for (int i = 0; i < count; i++) a *= function(i);
            return a;
        }
        #endregion

        #region complex[]
        public static complex[] Average(int count, Func<int, complex[]> function) { return Sum(count, function).LetDiv(count); }
        public static complex[] Sum(int count, Func<int, complex[]> function)
        {
            var a = (count == 0) ? default(complex[]) : function(0).CloneX();
            for (int i = 1; i < count; i++) a.LetAdd(function(i));
            return a;
        }
        #endregion

        #region complex[,]
        public static complex[,] Average(int count, Func<int, complex[,]> function) { return Sum(count, function).LetDiv(count); }
        public static complex[,] Sum(int count, Func<int, complex[,]> function)
        {
            var a = (count == 0) ? default(complex[,]) : function(0).CloneX();
            for (int i = 1; i < count; i++) a.LetAdd(function(i));
            return a;
        }
        #endregion

        #region Variance
        public static double StandardErrorMean<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).StandardErrorMean(); }
        public static double StandardErrorMean(this IEnumerable<double> source)
        {
            var r = source._Variance();
            if (r.v0 < 2) ThrowException.InvalidOperation("count < 2");
            return Math.Sqrt(r.v1 / (r.v0 - 1) / r.v0);
        }
        public static double StandardDeviation<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).StandardDeviation(); }
        public static double StandardDeviation(this IEnumerable<double> source) { return Math.Sqrt(source.Variance()); }
        public static double Variance<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).Variance(); }
        public static double Variance(this IEnumerable<double> source)
        {
            var r = source._Variance();
            if (r.v0 < 2) ThrowException.InvalidOperation("count < 2");
            return r.v1 / (r.v0 - 1);
        }
        public static double VariancePopulation<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).VariancePopulation(); }
        public static double VariancePopulation(this IEnumerable<double> source)
        {
            var r = source._Variance();
            if (r.v0 < 1) ThrowException.InvalidOperation("count < 1");
            return r.v1 / r.v0;
        }
        static V2<int, double> _Variance(this IEnumerable<double> source)
        {
            int count = 0;
            double moment1 = 0;
            double moment2 = 0;
            foreach (var element in source)
            {
                count++;
                checked
                {
                    moment1 += element;
                    moment2 += element * element;
                }
            }
            return New.V2(count, count == 0 ? 0 : Math.Max(0, moment2 - moment1 * moment1 / count));
        }
        #endregion

        #region Median
        public static double Median(this IList<double> data) { return Quantile(data, 0.5); }
        public static double MedianSorted(this IList<double> data) { return QuantileSorted(data, 0.5); }
        public static double Quantile(this IList<double> data, double quantile)
        {
            var index = data.SortIndex();
            return _QuantileSorted(index.Length, i => data[index[i]], quantile);
        }
        public static double QuantileSorted(this IList<double> data, double quantile)
        {
            return _QuantileSorted(data.Count, i => data[i], quantile);
        }
        static double _QuantileSorted(int count, Func<int, double> function, double quantile)
        {
            if (count == 0) ThrowException.Argument("count");
            if (function == null) ThrowException.Argument("function");
            if (quantile < 0 || quantile > 1) ThrowException.Argument("quantile");
            double r = (count - 1) * quantile;
            int index = (int)r;
            r -= index;
            return r == 0 ? function(index) : function(index) * (1 - r) + function(index + 1) * r;
        }
        #endregion

        #region Cumurative
        public static IEnumerable<int> Cumurative(this IEnumerable<int> source)
        {
            var a = default(int);
            return source.Select(element => checked(a += element));
        }
        public static IEnumerable<int> CumurativeZero(this IEnumerable<int> source)
        {
            var a = default(int);
            return source.Select(element => { var s = a; checked { a += element; } return s; });
        }
        public static IEnumerable<double> Cumurative(this IEnumerable<double> source)
        {
            var a = default(double);
            return source.Select(element => checked(a += element));
        }
        public static IEnumerable<double> CumurativeZero(this IEnumerable<double> source)
        {
            var a = default(double);
            return source.Select(element => { var s = a; checked { a += element; } return s; });
        }
        #endregion
        #endregion

        #region Integer functions
        static List<uint> PrimeNumbersUInt = new List<uint>() { 2, 3 };
        static uint PrimeNumbersUIntExamined = 3;
        static bool FindNextPrimeNumberUInt()
        {
            var v = PrimeNumbersUIntExamined;
            while (v != uint.MaxValue)
            {
                v += 2;
                var limit = (uint)Math.Sqrt(v);
                for (int i = 1; ; i++)
                {
                    var prime = PrimeNumbersUInt[i];
                    if (prime <= limit)
                    {
                        if (v % prime != 0) continue;
                        break;
                    }
                    PrimeNumbersUInt.Add(v);
                    PrimeNumbersUIntExamined = v;
                    return true;
                }
            }
            PrimeNumbersUIntExamined = v;
            return false;
        }

        static Tuple<int, int>[] FactorInteger_(int value)
        {
            lock (PrimeNumbersUInt)
            {
                var list = new List<Int2>();
                uint v = (uint)value;
                if (value == 0) { list.Add(new Int2(0, 1)); v = 1; }
                if (value < 0) { list.Add(new Int2(-1, 1)); v = (uint)-value; }
                for (int i = 0; v != 1; i++)
                {
                    if (PrimeNumbersUInt.Count <= i) FindNextPrimeNumberUInt();
                    var prime = PrimeNumbersUInt[i];
                    int count = 0;
                    while (v % prime == 0) { v /= prime; count++; }
                    if (count > 0) { list.Add(new Int2((int)prime, count)); continue; }
                    if (v < prime * prime) { list.Add(new Int2((int)v, 1)); break; }
                }
                return list.Select(p => Tuple.Create(p.v0, p.v1)).ToArray();
            }
        }
        public static Func<int, Tuple<int, int>[]> FactorInteger = New.Cache<int, Tuple<int, int>[]>(FactorInteger_);
        static int[] FactorIntegerExpanded_(int value)
        {
            return FactorInteger(value).SelectMany(p => Enumerable.Repeat(p.Item1, p.Item2)).ToArray();
        }
        public static Func<int, int[]> FactorIntegerExpanded = New.Cache<int, int[]>(FactorIntegerExpanded_);

        public static BigInteger MultinomialInteger(IEnumerable<int> source)
        {
            int total = 0;
            BigInteger product = 1;
            foreach (int element in source)
            {
                if (element < 0) ThrowException.Argument("element");
                total += element;
                product *= FactorialInteger(element);
            }
            return FactorialInteger(total) / product;
        }

        static List<BigInteger> FactorialIntegerBuffer = new List<BigInteger>() { 1 };
        static BigInteger FactorialInteger_(int value)
        {
            BigInteger product = FactorialIntegerBuffer.Last();
            for (int i = FactorialIntegerBuffer.Count; i <= value; i++)
            {
                product *= i;
                FactorialIntegerBuffer.Add(product);
            }
            return product;
        }
        public static BigInteger FactorialInteger(int value)
        {
            if (value < 0) ThrowException.ArgumentOutOfRange("value");
            return value < FactorialIntegerBuffer.Count ? FactorialIntegerBuffer[value] : FactorialInteger_(value);
        }

        public static BigInteger GreatestCommonDivisor(BigInteger value0, BigInteger value1)
        {
            while (true)
            {
                if (value0 < value1) Mt.Swap(ref value0, ref value1);
                var z = BigInteger.Remainder(value0, value1);
                if (z == 0) break;
                value0 = z;
            }
            return value1;
        }
        #endregion

        #region Special functions
        #region gamma functions
        readonly static double[] LogGammaCoefficients = { 57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199, .339946499848118887e-4, .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3, -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3, .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5 };
        public static double LogGamma(double value)
        {
            if (value <= 0) return double.NaN;
            double ser = 0.999999999999997092;
            for (int j = 0; j < LogGammaCoefficients.Length; j++) ser += LogGammaCoefficients[j] / (j + 1 + value);
            double tmp = value + 5.24218750000000000;
            return (value + 0.5) * Math.Log(tmp) - tmp + Math.Log(2.5066282746310005 * ser / value);
        }
        public static double Gamma(double value) { return Math.Exp(LogGamma(value)); }

        // 階乗
        // 22! まではdoubleで正確に表現可能
        // 170! まではdoubleで近似的に表現可能
        static double[] FactorialBuffer;
        static void Factorial_()
        {
            FactorialBuffer = new double[171];
            double f = 1;
            FactorialBuffer[0] = f;
            for (int i = 1; i < FactorialBuffer.Length; i++)
            {
                f *= i;
                FactorialBuffer[i] = f;
            }
        }
        public static double Factorial(int value)
        {
            if (FactorialBuffer == null) Factorial_();
            if (value < 0 || value >= FactorialBuffer.Length) ThrowException.ArgumentOutOfRange("value");
            return FactorialBuffer[value];
        }

        static double[] LogFactorialBuffer;
        static void LogFactorial_()
        {
            LogFactorialBuffer = new double[2000];
            double f = 0;
            for (int i = 2; i < LogFactorialBuffer.Length; i++)
            {
                f += Math.Log(i);
                LogFactorialBuffer[i] = i <= 22 ? Math.Log(Factorial(i)) : f;
            }
        }
        public static double LogFactorial(int value)
        {
            if (value < 0) ThrowException.ArgumentOutOfRange("value");
            if (LogFactorialBuffer == null) LogFactorial_();
            if (value < LogFactorialBuffer.Length) return LogFactorialBuffer[value];
            return Mt.LogGamma(value + 1.0);
        }

        // same to Pochhammer function
        public static double RisingFactorial(int value, int count)
        {
            if (count == 0) return 1;
            if (count == 1) return value;
            int s, v0, v1;
            if (value > 0)
            {
                s = 1;
                v0 = value + count - 1;
                v1 = value - 1;
            }
            else
            {
                s = 1 - 2 * (count & 1);
                v0 = -value;
                v1 = -value - count;
            }
            if (v0 < 0) return double.NaN;
            if (v1 < 0) return 0;
            return s * (v1 < 171 && v0 < 171 ? Factorial(v0) / Factorial(v1) : Math.Exp(LogFactorial(v0) - LogFactorial(v1)));
        }
        // FactorialPowerと同じ
        public static double FallingFactorial(int value, int count)
        {
            if (count == 0) return 1;
            if (count == 1) return value;
            int s, v0, v1;
            if (value >= count)
            {
                s = 1;
                v0 = value;
                v1 = value - count;
            }
            else
            {
                s = 1 - 2 * (count & 1);
                v0 = -value + count - 1;
                v1 = -value - 1;
            }
            if (v0 < 0) return double.NaN;
            if (v1 < 0) return 0;
            return s * (v1 < 171 && v0 < 171 ? Factorial(v0) / Factorial(v1) : Math.Exp(LogFactorial(v0) - LogFactorial(v1)));
        }

        public static double LogBinomial(int left, int right)
        {
            if (left < 0) ThrowException.ArgumentOutOfRange("left");
            int val2 = left - right;
            if (right < 0 || val2 < 0) return double.NegativeInfinity;
            if (right == 0 || val2 == 0) return 0.0;
            return Mt.LogFactorial(left) - (Mt.LogFactorial(right) + Mt.LogFactorial(val2));
        }
        public static double Binomial(int left, int right)
        {
            if (left < 0) ThrowException.ArgumentOutOfRange("left");
            int val2 = left - right;
            if (right < 0 || val2 < 0) return 0.0;
            if (right == 0 || val2 == 0) return 1.0;
            return Math.Round(Math.Exp(Mt.LogFactorial(left) - (Mt.LogFactorial(right) + Mt.LogFactorial(val2))), MidpointRounding.AwayFromZero);
        }

        public static double LogMultinomial(IEnumerable<int> source)
        {
            int total = 0;
            double sum = 0;
            foreach (int element in source)
            {
                if (element < 0) ThrowException.ArgumentOutOfRange("element");
                total += element;
                sum += Mt.LogFactorial(element);
            }
            return Mt.LogFactorial(total) - sum;
        }
        public static double Multinomial(IEnumerable<int> table)
        {
            return Math.Round(Math.Exp(LogMultinomial(table)), MidpointRounding.AwayFromZero);
        }

        // beta function
        public static double LogBeta(double value0, double value1)
        {
            return Mt.LogGamma(value0) + Mt.LogGamma(value1) - Mt.LogGamma(value0 + value1);
        }
        public static double Beta(double value0, double value1)
        {
            return Math.Exp(LogBeta(value0, value1));
        }
        #endregion

        #region incomplete gamma functions
        readonly static double[] Gauleg18y = new double[] {
            0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
            0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
            0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483,
            0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
            0.87126389619061517, 0.95698180152629142
        };
        readonly static double[] Gauleg18w = new double[] {
            0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
            0.034213810770299537, 0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
            0.058860144245324798, 0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
            0.076598410645870640, 0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
            0.085346685739338721, 0.085983275670394821
        };
        static double gammpapprox(double value, double gamma, bool psig)
        {
            double a1 = gamma - 1, lna1 = Math.Log(a1), sqrta1 = Math.Sqrt(a1);
            double xu = (value > a1) ?
                Math.Max(a1 + 11.5 * sqrta1, value + 6 * sqrta1) :
                Math.Max(0, Math.Min(a1 - 7.5 * sqrta1, value - 5 * sqrta1));
            double sum = 0;
            for (int j = 0; j < Gauleg18y.Length; j++)
            {
                double t = value + (xu - value) * Gauleg18y[j];
                sum += Gauleg18w[j] * Math.Exp(a1 - t + a1 * (Math.Log(t) - lna1));
            }
            double ans = sum * (xu - value) * Math.Exp(a1 * (lna1 - 1) - Mt.LogGamma(gamma));
            return psig ? (ans > 0 ? 1 : 0) - ans : (ans >= 0 ? 0 : 1) + ans;
        }
        static double gser(double value, double gamma)
        {
            double sum = 1 / gamma;
            double del = sum;
            for (int i = 1; ; i++)
            {
                del *= value / (gamma + i);
                sum += del;
                if (Math.Abs(del) < Math.Abs(sum) * Mt.DoubleEps) break;
            }
            return sum * Math.Exp(gamma * Math.Log(value) - value - Mt.LogGamma(gamma));
        }
        static double gcf(double value, double gamma)
        {
            double c = double.PositiveInfinity;
            double d = 1 / (value - gamma + 1);
            double h = d;
            for (int n = 1; ; n++)
            {
                double an = n * (gamma - n);
                double bn = (value - gamma + 1) + n * 2;
                c = bn + an / c; if (Math.Abs(c) < DoubleFpMin) c = DoubleFpMin;
                d = bn + an * d; if (Math.Abs(d) < DoubleFpMin) d = DoubleFpMin;
                if (c == d) break;
                h *= c / d;
                d = 1.0 / d;
            }
            return h * Math.Exp(gamma * Math.Log(value) - value - Mt.LogGamma(gamma));
        }
        // integrate [0, value] t^(gamma-1) e^-t dt / Gamma(gamma)
        public static double IncompleteGammaLower(double value, double gamma)
        {
            if (gamma <= 0) return double.NaN;
            if (value == double.PositiveInfinity) return 1;
            if (value <= 0) return 0;
            if (gamma >= 100) return gammpapprox(value, gamma, true);
            return (value < gamma + 1) ? gser(value, gamma) : 1 - gcf(value, gamma);
        }
        // integrate [value, inf) t^(gamma-1) e^-t dt / Gamma(gamma)
        public static double IncompleteGammaUpper(double value, double gamma)
        {
            if (gamma <= 0) return double.NaN;
            if (value == double.PositiveInfinity) return 0;
            if (value <= 0) return 1;
            if (gamma >= 100.0) return gammpapprox(value, gamma, false);
            return (value < gamma + 1) ? 1 - gser(value, gamma) : gcf(value, gamma);
        }
        static double InverseIncompleteGamma(double value, double gamma)
        {
            if (gamma <= 0) return double.NaN;
            if (value >= 1) return Math.Max(100, gamma + 100 * Math.Sqrt(gamma));
            if (value <= 0) return 0;

            const double eps = 1e-8;
            double x, t, lna1 = 0, afac = 0, a1 = gamma - 1;
            double gln = Mt.LogGamma(gamma);
            if (gamma > 1)
            {
                lna1 = Math.Log(a1);
                afac = Math.Exp(a1 * (lna1 - 1) - gln);
                t = Math.Sqrt(-2 * Math.Log((value < 0.5) ? value : 1 - value));
                x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;
                if (value < 0.5) x = -x;
                x = Math.Max(1.0e-3, gamma * Math.Pow(1 - 1 / (9 * gamma) - x / (3 * Math.Sqrt(gamma)), 3));
            }
            else
            {
                t = 1 - gamma * (0.253 + gamma * 0.12);
                x = (value < t) ?
                    Math.Pow(value / t, 1 / gamma) :
                    1 - Math.Log(1 - (value - t) / (1 - t));
            }
            for (int j = 0; j < 12; j++)
            {
                if (x <= 0) return 0;
                double err = Mt.IncompleteGammaLower(x, gamma) - value;
                t = (gamma > 1) ?
                    afac * Math.Exp(-(x - a1) + a1 * (Math.Log(x) - lna1)) :
                    Math.Exp(-x + a1 * Math.Log(x) - gln);
                double u = err / t;
                x -= (t = u / (1 - 0.5 * Math.Min(1, u * (a1 / x - 1))));
                if (x <= 0) x = 0.5 * (x + t);
                if (Math.Abs(t) < eps * x) break;
            }
            return x;
        }
        public static double Erf(double value)
        {
            return Mt.IncompleteGammaLower(value * value, 0.5) * Math.Sign(value);
        }

        // normal distribution
        static public double StandardNormalDistribution(double value)
        {
            return Math.Exp(value * value * -0.5) * Sqrt2PIinv;
        }
        static public double StandardNormalDistributionLower(double value)
        {
            return (1 + Mt.Erf(value * Sqrt2inv)) / 2;
        }
        static public double StandardNormalDistributionUpper(double value)
        {
            return (1 - Mt.Erf(value * Sqrt2inv)) / 2;
        }
        static public double NormalDistribution(double value, double mean, double variance)
        {
            if (variance <= 0) ThrowException.Argument("variance");
            return Math.Exp(Mt.Sq(value - mean) / variance * -0.5) / Math.Sqrt(Mt.PI2 * variance);
        }

        // chi-square distribution
        public static double ChiSquareDistribution(double value, double freedom)
        {
            if (freedom < 0) ThrowException.Argument("freedom");
            if (value < 0) return 0;
            if (value == 0)
            {
                if (freedom < 2) return double.PositiveInfinity;
                if (freedom == 2) return 0.5;
                return 0;
            }
            double f = freedom / 2;
            return Math.Exp(-0.5 * value + (f - 1) * Math.Log(value) - f * Mt.Ln2 - Mt.LogGamma(f));
        }
        public static double ChiSquareDistributionLower(double value, double freedom)
        {
            return Mt.IncompleteGammaLower(value / 2, freedom / 2);
        }
        public static double ChiSquareDistributionUpper(double value, double freedom)
        {
            return Mt.IncompleteGammaUpper(value / 2, freedom / 2);
        }
        #endregion

        #region incomplete beta functions
        static double AbsMaxFPMIN(double value)
        {
            return Math.Abs(value) < DoubleFpMin ? DoubleFpMin : value;
        }
        static double betacf(double a, double b, double x)
        {
            double qab = a + b;
            double qap = a + 1;
            double qam = a - 1;
            double c = 1;
            double d = 1 / AbsMaxFPMIN(1 - qab * x / qap);
            double h = d;
            for (int m = 1; m < 10000; m++)
            {
                int m2 = 2 * m;
                double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1 / AbsMaxFPMIN(1 + aa * d);
                c = AbsMaxFPMIN(1 + aa / c);
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = 1 / AbsMaxFPMIN(1 + aa * d);
                c = AbsMaxFPMIN(1 + aa / c);
                h *= d * c;
                if (Math.Abs(d * c - 1) <= DoubleEps) break;
            }
            return h;
        }

        static double betaiapprox(double a, double b, double x)
        {
            double xu;
            double mu = a / (a + b);
            double lnmu = Math.Log(mu);
            double lnmuc = Math.Log(1 - mu);
            double t = Math.Sqrt(a * b / (Mt.Sq(a + b) * (a + b + 1)));
            if (x > a / (a + b))
            {
                if (x >= 1) return 1;
                xu = Math.Min(1.0, Math.Max(mu + 10 * t, x + 5 * t));
            }
            else
            {
                if (x <= 0) return 0;
                xu = Math.Max(0.0, Math.Min(mu - 10 * t, x - 5 * t));
            }
            double sum = 0;
            for (int j = 0; j < 18; j++)
            {
                t = x + (xu - x) * Gauleg18y[j];
                sum += Gauleg18w[j] * Math.Exp((a - 1) * (Math.Log(t) - lnmu) + (b - 1) * (Math.Log(1 - t) - lnmuc));
            }
            double ans = sum * (xu - x) * Math.Exp((a - 1) * lnmu - Mt.LogGamma(a) + (b - 1) * lnmuc - Mt.LogGamma(b) + Mt.LogGamma(a + b));
            return ans > 0 ? 1 - ans : -ans;
        }
        public static double IncompleteBeta(double value, double param1, double param2)
        {
            if (value < 0 || value > 1) ThrowException.ArgumentOutOfRange("value");
            if (param1 <= 0) ThrowException.ArgumentOutOfRange("param1");
            if (param2 <= 0) ThrowException.ArgumentOutOfRange("param2");
            if (value == 0 || value == 1) return value;
            if (param1 > 3000 && param2 > 3000) return betaiapprox(param1, param2, value);
            double bt = Math.Exp(
                Mt.LogGamma(param1 + param2) - Mt.LogGamma(param1) - Mt.LogGamma(param2)
                + param1 * Math.Log(value) + param2 * Math.Log(1 - value));
            return (value < (param1 + 1) / (param1 + param2 + 2)) ?
                bt * betacf(param1, param2, value) / param1 :
                1 - bt * betacf(param2, param1, 1 - value) / param2;
        }

        // Student t distribution
        static public double StudentTDistribution(double value, double freedom)
        {
            double n2 = (freedom + 1) / 2;
            return Mt.Gamma(n2) / (Mt.Gamma(freedom / 2) * Math.Sqrt(Math.PI * freedom)) * Math.Pow(1 + value * value / freedom, -n2);
        }
        // ∫∞~-t and t~∞
        static public double StudentTDistributionBilateral(double value, double freedom)
        {
            return Mt.IncompleteBeta(freedom / (freedom + value * value), freedom / 2, 0.5);
        }
        // ∫t~∞
        static public double StudentTDistributionUpper(double value, double freedom)
        {
            double a = Mt.StudentTDistributionBilateral(value, freedom) / 2;
            return value < 0 ? 1 - a : a;
        }
        // ∫-∞~t
        static public double StudentTDistributionLower(double value, double freedom)
        {
            double a = Mt.StudentTDistributionBilateral(value, freedom) / 2;
            return value < 0 ? a : 1 - a;
        }

        // F distribution
        // ∫f~∞, Q(F|v1,v2) <0.01 で有意
        static public double FDistributionUpper(double value, double freedom1, double freedom2)
        {
            return Mt.IncompleteBeta(freedom2 / (freedom2 + freedom1 * value), freedom2 / 2, freedom1 / 2);
        }
        // ∫-∞~f = ∫0~f
        static public double FDistributionLower(double value, double freedom1, double freedom2)
        {
            return 1 - Mt.FDistributionUpper(value, freedom1, freedom2);
        }
        #endregion

        #region elliptic functions
        public static double EllipticTheta3(double phase, double radius)
        {
            if (radius < 0 || radius >= 1) ThrowException.ArgumentOutOfRange("radius");
            if (radius == 0) return 1;
            int digit = (int)Math.Ceiling(Math.Sqrt(Math.Log(DoubleEps) / Math.Log(radius)));
            double a = 0;
            for (int d = 1; d <= digit; d++)
                a += Math.Cos(2 * d * phase) * Math.Pow(radius, d * d);
            return 1 + 2 * a;
        }
        #endregion
        #endregion

        #region Linear functions
        #region fixed size
        public static double Norm1(this Double2 array) { return Math.Abs(array.X) + Math.Abs(array.Y); }
        public static double Norm1(this Double3 array) { return Math.Abs(array.X) + Math.Abs(array.Y) + Math.Abs(array.Z); }
        public static double SqNorm2(this Double2 array) { return array.X * array.X + array.Y * array.Y; }
        public static double SqNorm2(this Double3 array) { return array.X * array.X + array.Y * array.Y + array.Z * array.Z; }
        public static double Norm2(this Double2 array) { return Math.Sqrt(array.SqNorm2()); }
        public static double Norm2(this Double3 array) { return Math.Sqrt(array.SqNorm2()); }
        public static double Inner(this Double2 array0, Double2 array1) { return array0.X * array1.X + array0.Y * array1.Y; }
        public static double Inner(this Double3 array0, Double3 array1) { return array0.X * array1.X + array0.Y * array1.Y + array0.Z * array1.Z; }
        public static Double2 NormalizeNorm2(Double2 array) { return array / Norm2(array); }
        public static Double3 NormalizeNorm2(Double3 array) { return array / Norm2(array); }
        public static Double3 Outer(Double3 array0, Double3 array1)
        {
            return new Double3(
                array0.Y * array1.Z - array0.Z * array1.Y,
                array0.Z * array1.X - array0.X * array1.Z,
                array0.X * array1.Y - array0.Y * array1.X
            );
        }
        #endregion

        #region double[]
        static TSource0[] ToZero<TSource0, TSource1>(TSource0[] array0, TSource1[] array1) { Ex.SizeCheck(array0, array1); return array0.ZeroTo(); }

        public unsafe static double Sum(this double[] array, Func<double, double> selector) { fixed (double* p = array) return Us.Sum(p, array.Length, selector); }
        public unsafe static double Sum(this double[] array) { fixed (double* p = array) return Us.Sum(p, array.Length); }
        public unsafe static double Average(this double[] array, Func<double, double> selector) { return Sum(array, selector) / array.Length; }
        public unsafe static double Average(this double[] array) { return Sum(array) / array.Length; }
        public unsafe static double Norm(this double[] array, double nu) { fixed (double* p = array) return Math.Pow(Us.SumPow(p, array.Length, nu), 1 / nu); }
        public unsafe static double Norm1(this double[] array) { fixed (double* p = array) return Us.SumAbs(p, array.Length); }
        public unsafe static double Norm2(this double[] array) { return Math.Sqrt(SqNorm2(array)); }
        public unsafe static double SqNorm2(this double[] array) { fixed (double* p = array) return Us.SumSq(p, array.Length); }
        public unsafe static double Min(this double[] array) { fixed (double* p = array) return Us.Min(p, array.Length); }
        public unsafe static double Max(this double[] array) { fixed (double* p = array) return Us.Max(p, array.Length); }
        public unsafe static double MaxAbs(this double[] array) { fixed (double* p = array) return Us.MaxAbs(p, array.Length); }
        public unsafe static double MaxSq(this double[] array) { fixed (double* p = array) return Us.MaxSq(p, array.Length); }
        public unsafe static double Inner(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static double NormSub(double[] array0, double[] array1, double nu) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumPowSub(p, q, array0.Length, nu); }
        public unsafe static double Norm1Sub(double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumAbsSub(p, q, array0.Length); }
        public unsafe static double Norm2Sub(double[] array0, double[] array1) { return Math.Sqrt(SqNorm2Sub(array0, array1)); }
        public unsafe static double SqNorm2Sub(double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumSqSub(p, q, array0.Length); }

        public unsafe static double[] LetNeg(this double[] array) { fixed (double* p = array) Us.Neg(p, p, array.Length); return array; }
        public unsafe static double[] LetRev(this double[] array) { fixed (double* p = array) Us.LetRev(p, array.Length); return array; }
        public unsafe static double[] Let(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static double[] LetAdd(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static double[] LetSub(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static double[] LetMul(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static double[] LetDiv(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static double[] LetMod(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mod(p, p, q, array0.Length); return array0; }
        public unsafe static double[] LetSubr(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static double[] LetDivr(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static double[] LetModr(this double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mod(p, q, p, array0.Length); return array0; }
        public unsafe static double[] LetAdd(this double[] array, double scalar) { if (scalar != 0) fixed (double* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static double[] LetSub(this double[] array, double scalar) { if (scalar != 0) fixed (double* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static double[] LetMul(this double[] array, double scalar) { if (scalar != 1) fixed (double* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static double[] LetDiv(this double[] array, double scalar) { if (scalar != 1) fixed (double* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static double[] LetMod(this double[] array, double scalar) { fixed (double* p = array) Us.Mod(p, p, scalar, array.Length); return array; }
        public unsafe static double[] LetSubr(this double[] array, double scalar) { fixed (double* p = array) Us.Sub(p, scalar, p, array.Length); return array; }
        public unsafe static double[] LetDivr(this double[] array, double scalar) { fixed (double* p = array) Us.Div(p, scalar, p, array.Length); return array; }
        public unsafe static double[] LetModr(this double[] array, double scalar) { fixed (double* p = array) Us.Mod(p, scalar, p, array.Length); return array; }
        public unsafe static double[] Pos(this double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static double[] Neg(this double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Neg(r, p, result.Length); return result; }
        public unsafe static double[] Rev(this double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Rev(r, p, result.Length); return result; }
        public unsafe static double[] Add(this double[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static double[] Sub(this double[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static double[] Mul(this double[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static double[] Div(this double[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static double[] Mod(this double[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Mod(r, p, q, result.Length); return result; }
        public unsafe static double[] Add(this double[] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static double[] Sub(this double[] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static double[] Mul(this double[] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static double[] Div(this double[] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static double[] Mod(this double[] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mod(r, p, scalar, result.Length); return result; }
        public unsafe static double[] Add(double scalar, double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static double[] Sub(double scalar, double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static double[] Mul(double scalar, double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static double[] Div(double scalar, double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }
        public unsafe static double[] Mod(double scalar, double[] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mod(r, scalar, p, result.Length); return result; }

        public unsafe static double[] LetAddMul(this double[] array0, double[] array1, double scalar) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.LetAddMul(p, q, scalar, array0.Length); return array0; }
        public unsafe static double[] LetMulAdd(this double[] array0, double scalar, double[] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.LetMulAdd(p, scalar, q, array0.Length); return array0; }
        public unsafe static double[] AddMul(double[] array0, double[] array1, double scalar) { Ex.SizeCheck(array0, array1); var result = array0.CloneX(); fixed (double* r = result, p = array1) Us.LetAddMul(r, p, scalar, result.Length); return result; }
        public unsafe static double[] LetNormalizeSum(this double[] data) { data.LetDiv(data.Sum()); return data; }
        public unsafe static double[] LetNormalizeNorm1(this double[] data) { data.LetDiv(data.Norm1()); return data; }
        public unsafe static double[] LetNormalizeNorm2(this double[] data) { data.LetDiv(data.Norm2()); return data; }
        #endregion

        #region double[,]
        static TSource0[,] ToZero<TSource0, TSource1>(TSource0[,] array0, TSource1[,] array1) { Ex.SizeCheck(array0, array1); return array0.ZeroTo(); }

        public unsafe static double Sum(this double[,] array, Func<double, double> selector) { fixed (double* p = array) return Us.Sum(p, array.Length, selector); }
        public unsafe static double Sum(this double[,] array) { fixed (double* p = array) return Us.Sum(p, array.Length); }
        public unsafe static double Average(this double[,] array, Func<double, double> selector) { return Sum(array, selector) / array.Length; }
        public unsafe static double Average(this double[,] array) { return Sum(array) / array.Length; }
        public unsafe static double Norm(this double[,] array, double nu) { fixed (double* p = array) return Math.Pow(Us.SumPow(p, array.Length, nu), 1 / nu); }
        public unsafe static double Norm1(this double[,] array) { fixed (double* p = array) return Us.SumAbs(p, array.Length); }
        public unsafe static double Norm2(this double[,] array) { return Math.Sqrt(SqNorm2(array)); }
        public unsafe static double SqNorm2(this double[,] array) { fixed (double* p = array) return Us.SumSq(p, array.Length); }
        public unsafe static double Min(this double[,] array) { fixed (double* p = array) return Us.Min(p, array.Length); }
        public unsafe static double Max(this double[,] array) { fixed (double* p = array) return Us.Max(p, array.Length); }
        public unsafe static double MaxAbs(this double[,] array) { fixed (double* p = array) return Us.MaxAbs(p, array.Length); }
        public unsafe static double MaxSq(this double[,] array) { fixed (double* p = array) return Us.MaxSq(p, array.Length); }
        public unsafe static double Inner(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static double NormSub(double[,] array0, double[,] array1, double nu) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumPowSub(p, q, array0.Length, nu); }
        public unsafe static double Norm1Sub(double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumAbsSub(p, q, array0.Length); }
        public unsafe static double Norm2Sub(double[,] array0, double[,] array1) { return Math.Sqrt(SqNorm2Sub(array0, array1)); }
        public unsafe static double SqNorm2Sub(double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumSqSub(p, q, array0.Length); }

        public unsafe static double[,] LetNeg(this double[,] array) { fixed (double* p = array) Us.Neg(p, p, array.Length); return array; }
        public unsafe static double[,] LetRev(this double[,] array) { fixed (double* p = array) Us.LetRev(p, array.Length); return array; }
        public unsafe static double[,] Let(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static double[,] LetAdd(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static double[,] LetSub(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static double[,] LetMul(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static double[,] LetDiv(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static double[,] LetMod(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mod(p, p, q, array0.Length); return array0; }
        public unsafe static double[,] LetSubr(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static double[,] LetDivr(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static double[,] LetModr(this double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mod(p, q, p, array0.Length); return array0; }
        public unsafe static double[,] LetAdd(this double[,] array, double scalar) { if (scalar != 0) fixed (double* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static double[,] LetSub(this double[,] array, double scalar) { if (scalar != 0) fixed (double* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static double[,] LetMul(this double[,] array, double scalar) { if (scalar != 1) fixed (double* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static double[,] LetDiv(this double[,] array, double scalar) { if (scalar != 1) fixed (double* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static double[,] LetMod(this double[,] array, double scalar) { fixed (double* p = array) Us.Mod(p, p, scalar, array.Length); return array; }
        public unsafe static double[,] LetSubr(this double[,] array, double scalar) { fixed (double* p = array) Us.Sub(p, scalar, p, array.Length); return array; }
        public unsafe static double[,] LetDivr(this double[,] array, double scalar) { fixed (double* p = array) Us.Div(p, scalar, p, array.Length); return array; }
        public unsafe static double[,] LetModr(this double[,] array, double scalar) { fixed (double* p = array) Us.Mod(p, scalar, p, array.Length); return array; }
        public unsafe static double[,] Pos(this double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static double[,] Neg(this double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Neg(r, p, result.Length); return result; }
        public unsafe static double[,] Rev(this double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Rev(r, p, result.Length); return result; }
        public unsafe static double[,] Add(this double[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static double[,] Sub(this double[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static double[,] Mul(this double[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static double[,] Div(this double[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static double[,] Mod(this double[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Mod(r, p, q, result.Length); return result; }
        public unsafe static double[,] Add(this double[,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static double[,] Sub(this double[,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static double[,] Mul(this double[,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static double[,] Div(this double[,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static double[,] Mod(this double[,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mod(r, p, scalar, result.Length); return result; }
        public unsafe static double[,] Add(double scalar, double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static double[,] Sub(double scalar, double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static double[,] Mul(double scalar, double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static double[,] Div(double scalar, double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }
        public unsafe static double[,] Mod(double scalar, double[,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mod(r, scalar, p, result.Length); return result; }

        public unsafe static double[,] LetAddMul(this double[,] array0, double[,] array1, double scalar) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.LetAddMul(p, q, scalar, array0.Length); return array0; }
        public unsafe static double[,] LetMulAdd(this double[,] array0, double scalar, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.LetMulAdd(p, scalar, q, array0.Length); return array0; }
        public unsafe static double[,] AddMul(double[,] array0, double[,] array1, double scalar) { Ex.SizeCheck(array0, array1); var result = array0.CloneX(); fixed (double* r = result, p = array1) Us.LetAddMul(r, p, scalar, result.Length); return result; }
        #endregion

        #region double[, ,]
        static TSource0[, ,] ToZero<TSource0, TSource1>(TSource0[, ,] array0, TSource1[, ,] array1) { Ex.SizeCheck(array0, array1); return array0.ZeroTo(); }

        public unsafe static double Sum(this double[, ,] array, Func<double, double> selector) { fixed (double* p = array) return Us.Sum(p, array.Length, selector); }
        public unsafe static double Sum(this double[, ,] array) { fixed (double* p = array) return Us.Sum(p, array.Length); }
        public unsafe static double Average(this double[, ,] array, Func<double, double> selector) { return Sum(array, selector) / array.Length; }
        public unsafe static double Average(this double[, ,] array) { return Sum(array) / array.Length; }
        public unsafe static double Norm(this double[, ,] array, double nu) { fixed (double* p = array) return Math.Pow(Us.SumPow(p, array.Length, nu), 1 / nu); }
        public unsafe static double Norm1(this double[, ,] array) { fixed (double* p = array) return Us.SumAbs(p, array.Length); }
        public unsafe static double Norm2(this double[, ,] array) { return Math.Sqrt(SqNorm2(array)); }
        public unsafe static double SqNorm2(this double[, ,] array) { fixed (double* p = array) return Us.SumSq(p, array.Length); }
        public unsafe static double Min(this double[, ,] array) { fixed (double* p = array) return Us.Min(p, array.Length); }
        public unsafe static double Max(this double[, ,] array) { fixed (double* p = array) return Us.Max(p, array.Length); }
        public unsafe static double MaxAbs(this double[, ,] array) { fixed (double* p = array) return Us.MaxAbs(p, array.Length); }
        public unsafe static double MaxSq(this double[, ,] array) { fixed (double* p = array) return Us.MaxSq(p, array.Length); }
        public unsafe static double Inner(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static double NormSub(double[, ,] array0, double[, ,] array1, double nu) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumPowSub(p, q, array0.Length, nu); }
        public unsafe static double Norm1Sub(double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumAbsSub(p, q, array0.Length); }
        public unsafe static double Norm2Sub(double[, ,] array0, double[, ,] array1) { return Math.Sqrt(SqNorm2Sub(array0, array1)); }
        public unsafe static double SqNorm2Sub(double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) return Us.SumSqSub(p, q, array0.Length); }

        public unsafe static double[, ,] LetNeg(this double[, ,] array) { fixed (double* p = array) Us.Neg(p, p, array.Length); return array; }
        public unsafe static double[, ,] LetRev(this double[, ,] array) { fixed (double* p = array) Us.LetRev(p, array.Length); return array; }
        public unsafe static double[, ,] Let(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static double[, ,] LetAdd(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static double[, ,] LetSub(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static double[, ,] LetMul(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static double[, ,] LetDiv(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static double[, ,] LetMod(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mod(p, p, q, array0.Length); return array0; }
        public unsafe static double[, ,] LetSubr(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static double[, ,] LetDivr(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static double[, ,] LetModr(this double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.Mod(p, q, p, array0.Length); return array0; }
        public unsafe static double[, ,] LetAdd(this double[, ,] array, double scalar) { if (scalar != 0) fixed (double* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static double[, ,] LetSub(this double[, ,] array, double scalar) { if (scalar != 0) fixed (double* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static double[, ,] LetMul(this double[, ,] array, double scalar) { if (scalar != 1) fixed (double* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static double[, ,] LetDiv(this double[, ,] array, double scalar) { if (scalar != 1) fixed (double* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static double[, ,] LetMod(this double[, ,] array, double scalar) { fixed (double* p = array) Us.Mod(p, p, scalar, array.Length); return array; }
        public unsafe static double[, ,] LetSubr(this double[, ,] array, double scalar) { fixed (double* p = array) Us.Sub(p, scalar, p, array.Length); return array; }
        public unsafe static double[, ,] LetDivr(this double[, ,] array, double scalar) { fixed (double* p = array) Us.Div(p, scalar, p, array.Length); return array; }
        public unsafe static double[, ,] LetModr(this double[, ,] array, double scalar) { fixed (double* p = array) Us.Mod(p, scalar, p, array.Length); return array; }
        public unsafe static double[, ,] Pos(this double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static double[, ,] Neg(this double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Neg(r, p, result.Length); return result; }
        public unsafe static double[, ,] Rev(this double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Rev(r, p, result.Length); return result; }
        public unsafe static double[, ,] Add(this double[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static double[, ,] Sub(this double[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static double[, ,] Mul(this double[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static double[, ,] Div(this double[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static double[, ,] Mod(this double[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (double* r = result, p = array0, q = array1) Us.Mod(r, p, q, result.Length); return result; }
        public unsafe static double[, ,] Add(this double[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static double[, ,] Sub(this double[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static double[, ,] Mul(this double[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static double[, ,] Div(this double[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static double[, ,] Mod(this double[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mod(r, p, scalar, result.Length); return result; }
        public unsafe static double[, ,] Add(double scalar, double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static double[, ,] Sub(double scalar, double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static double[, ,] Mul(double scalar, double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static double[, ,] Div(double scalar, double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }
        public unsafe static double[, ,] Mod(double scalar, double[, ,] array) { var result = array.ZeroTo(); fixed (double* r = result, p = array) Us.Mod(r, scalar, p, result.Length); return result; }

        public unsafe static double[, ,] LetAddMul(this double[, ,] array0, double[, ,] array1, double scalar) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.LetAddMul(p, q, scalar, array0.Length); return array0; }
        public unsafe static double[, ,] LetMulAdd(this double[, ,] array0, double scalar, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (double* p = array0, q = array1) Us.LetMulAdd(p, scalar, q, array0.Length); return array0; }
        public unsafe static double[, ,] AddMul(double[, ,] array0, double[, ,] array1, double scalar) { Ex.SizeCheck(array0, array1); var result = array0.CloneX(); fixed (double* r = result, p = array1) Us.LetAddMul(r, p, scalar, result.Length); return result; }
        #endregion

        #region complex[]
        public unsafe static complex Sum(this complex[] array) { fixed (complex* p = array) return Us.Sum(p, array.Length); }
        public unsafe static complex Average(this complex[] array) { return array.Sum() / array.Length; }
        public unsafe static complex Inner(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static complex InnerConj(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumCul(p, q, array0.Length); }
        public unsafe static complex Inner(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static complex InnerConj(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) return Us.SumCul(p, q, array0.Length); }
        public unsafe static double Norm(this complex[] array, double nu) { fixed (complex* p = array) return Math.Pow(Us.SumPowAbs(p, array.Length, nu), 1 / nu); }
        public unsafe static double Norm1(this complex[] array) { fixed (complex* p = array) return Us.SumAbs(p, array.Length); }
        public unsafe static double Norm2(this complex[] array) { return Math.Sqrt(SqNorm2(array)); }
        public unsafe static double SqNorm2(this complex[] array) { fixed (complex* p = array) return Us.SumSqAbs(p, array.Length); }
        public unsafe static double MaxAbs(this complex[] array) { return Math.Sqrt(MaxSqAbs(array)); }
        public unsafe static double MaxSqAbs(this complex[] array) { fixed (complex* p = array) return Us.MaxSqAbs(p, array.Length); }
        public unsafe static double NormSub(complex[] array0, complex[] array1, double nu) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumPowSub(p, q, array0.Length, nu); }
        public unsafe static double Norm1Sub(complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumAbsSub(p, q, array0.Length); }
        public unsafe static double Norm2Sub(complex[] array0, complex[] array1) { return Math.Sqrt(SqNorm2Sub(array0, array1)); }
        public unsafe static double SqNorm2Sub(complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumSqSub(p, q, array0.Length); }

        public unsafe static complex[] LetConj(this complex[] array) { fixed (complex* p = array) Us.Conj(p, p, array.Length); return array; }
        public unsafe static complex[] LetNeg(this complex[] array) { fixed (complex* p = array) Us.Neg(p, p, array.Length); return array; }
        public unsafe static complex[] LetRev(this complex[] array) { fixed (complex* p = array) Us.LetRev(p, array.Length); return array; }
        public unsafe static complex[] Let(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static complex[] LetAdd(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetSub(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetMul(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetDiv(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetSubr(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static complex[] LetDivr(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static complex[] LetConjMul(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Cul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetMulConj(this complex[] array0, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Cul(p, q, p, array0.Length); return array0; }
        public unsafe static complex[] Let(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static complex[] LetAdd(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetSub(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetMul(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetDiv(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetConjMul(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Cul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[] LetSubr(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static complex[] LetDivr(this complex[] array0, double[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static complex[] LetAdd(this complex[] array, complex scalar) { if (scalar != 0) fixed (complex* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetSub(this complex[] array, complex scalar) { if (scalar != 0) fixed (complex* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetMul(this complex[] array, complex scalar) { if (scalar != 1) fixed (complex* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetDiv(this complex[] array, complex scalar) { if (scalar != 1) fixed (complex* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetConjMul(this complex[] array, complex scalar) { fixed (complex* p = array) Us.Cul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetSubr(this complex[] array, complex scalar) { fixed (complex* p = array) Us.Sub(p, scalar, p, array.Length); return array; }
        public unsafe static complex[] LetDivr(this complex[] array, complex scalar) { fixed (complex* p = array) Us.Div(p, scalar, p, array.Length); return array; }
        public unsafe static complex[] LetAdd(this complex[] array, double scalar) { if (scalar != 0) fixed (complex* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetSub(this complex[] array, double scalar) { if (scalar != 0) fixed (complex* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetMul(this complex[] array, double scalar) { if (scalar != 1) fixed (complex* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetDiv(this complex[] array, double scalar) { if (scalar != 1) fixed (complex* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] LetConjMul(this complex[] array, double scalar) { fixed (complex* p = array) Us.Cul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[] Tocomplex(double[] array0, double[] array1) { Ex.SizeCheck(array0, array1); var result = new complex[array0.Length]; fixed (complex* r = result) fixed (double* p = array0, q = array1) Us.Tocomplex(r, p, q, result.Length); return result; }
        public unsafe static complex[] Tocomplex(this double[] array) { var result = new complex[array.Length]; fixed (complex* r = result) fixed (double* p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static double[] Real(this complex[] array) { var result = new double[array.Length]; fixed (double* r = result) fixed (complex* p = array) Us.Real(r, p, result.Length); return result; }
        public unsafe static double[] Imag(this complex[] array) { var result = new double[array.Length]; fixed (double* r = result) fixed (complex* p = array) Us.Imag(r, p, result.Length); return result; }
        public unsafe static complex[] Conj(this complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Conj(r, p, result.Length); return result; }
        public unsafe static complex[] Pos(this complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static complex[] Neg(this complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Neg(r, p, result.Length); return result; }
        public unsafe static complex[] Rev(this complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Rev(r, p, result.Length); return result; }
        public unsafe static complex[] Add(this complex[] array0, complex[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[] Sub(this complex[] array0, complex[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static complex[] Mul(this complex[] array0, complex[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[] Div(this complex[] array0, complex[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static complex[] ConjMul(this complex[] array0, complex[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Cul(r, p, q, result.Length); return result; }
        public unsafe static complex[] Add(this complex[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[] Sub(this complex[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static complex[] Mul(this complex[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[] Div(this complex[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static complex[] ConjMul(this complex[] array0, double[] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Cul(r, p, q, result.Length); return result; }
        public unsafe static complex[] Add(this double[] array0, complex[] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[] Sub(this double[] array0, complex[] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Sub(r, q, p, result.Length); return result; }
        public unsafe static complex[] Mul(this double[] array0, complex[] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[] Div(this double[] array0, complex[] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Div(r, q, p, result.Length); return result; }
        public unsafe static complex[] Add(this complex[] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Sub(this complex[] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Mul(this complex[] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Div(this complex[] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] ConjMul(this complex[] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Cul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Add(complex scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Sub(complex scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static complex[] Mul(complex scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Div(complex scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }
        public unsafe static complex[] Add(this complex[] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Sub(this complex[] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Mul(this complex[] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Div(this complex[] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] ConjMul(this complex[] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Cul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Add(double scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Sub(double scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static complex[] Mul(double scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[] Div(double scalar, complex[] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }

        public unsafe static complex[] LetAddMul(this complex[] array0, complex[] array1, double scalar) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.LetAddMul(p, q, scalar, array0.Length); return array0; }
        public unsafe static complex[] LetMulAdd(this complex[] array0, double scalar, complex[] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.LetMulAdd(p, scalar, q, array0.Length); return array0; }
        public unsafe static complex[] AddMul(this complex[] array0, complex[] array1, double scalar) { Ex.SizeCheck(array0, array1); var result = array0.CloneX(); fixed (complex* r = result, p = array1) Us.LetAddMul(r, p, scalar, result.Length); return result; }
        #endregion

        #region complex[,]
        public unsafe static complex Sum(this complex[,] array) { fixed (complex* p = array) return Us.Sum(p, array.Length); }
        public unsafe static complex Average(this complex[,] array) { return array.Sum() / array.Length; }
        public unsafe static complex Inner(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static complex InnerConj(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumCul(p, q, array0.Length); }
        public unsafe static complex Inner(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static complex InnerConj(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) return Us.SumCul(p, q, array0.Length); }
        public unsafe static double Norm(this complex[,] array, double nu) { fixed (complex* p = array) return Math.Pow(Us.SumPowAbs(p, array.Length, nu), 1 / nu); }
        public unsafe static double Norm1(this complex[,] array) { fixed (complex* p = array) return Us.SumAbs(p, array.Length); }
        public unsafe static double Norm2(this complex[,] array) { return Math.Sqrt(SqNorm2(array)); }
        public unsafe static double SqNorm2(this complex[,] array) { fixed (complex* p = array) return Us.SumSqAbs(p, array.Length); }
        public unsafe static double MaxAbs(this complex[,] array) { return Math.Sqrt(MaxSqAbs(array)); }
        public unsafe static double MaxSqAbs(this complex[,] array) { fixed (complex* p = array) return Us.MaxSqAbs(p, array.Length); }
        public unsafe static double NormSub(complex[,] array0, complex[,] array1, double nu) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumPowSub(p, q, array0.Length, nu); }
        public unsafe static double Norm1Sub(complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumAbsSub(p, q, array0.Length); }
        public unsafe static double Norm2Sub(complex[,] array0, complex[,] array1) { return Math.Sqrt(SqNorm2Sub(array0, array1)); }
        public unsafe static double SqNorm2Sub(complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumSqSub(p, q, array0.Length); }

        public unsafe static complex[,] LetConj(this complex[,] array) { fixed (complex* p = array) Us.Conj(p, p, array.Length); return array; }
        public unsafe static complex[,] LetNeg(this complex[,] array) { fixed (complex* p = array) Us.Neg(p, p, array.Length); return array; }
        public unsafe static complex[,] LetRev(this complex[,] array) { fixed (complex* p = array) Us.LetRev(p, array.Length); return array; }
        public unsafe static complex[,] Let(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetAdd(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetSub(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetMul(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetDiv(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetSubr(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static complex[,] LetDivr(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static complex[,] LetConjMul(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Cul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetMulConj(this complex[,] array0, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Cul(p, q, p, array0.Length); return array0; }
        public unsafe static complex[,] Let(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetAdd(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetSub(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetMul(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetDiv(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetConjMul(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Cul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[,] LetSubr(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static complex[,] LetDivr(this complex[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static complex[,] LetAdd(this complex[,] array, complex scalar) { if (scalar != 0) fixed (complex* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetSub(this complex[,] array, complex scalar) { if (scalar != 0) fixed (complex* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetMul(this complex[,] array, complex scalar) { if (scalar != 1) fixed (complex* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetDiv(this complex[,] array, complex scalar) { if (scalar != 1) fixed (complex* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetConjMul(this complex[,] array, complex scalar) { fixed (complex* p = array) Us.Cul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetSubr(this complex[,] array, complex scalar) { fixed (complex* p = array) Us.Sub(p, scalar, p, array.Length); return array; }
        public unsafe static complex[,] LetDivr(this complex[,] array, complex scalar) { fixed (complex* p = array) Us.Div(p, scalar, p, array.Length); return array; }
        public unsafe static complex[,] LetAdd(this complex[,] array, double scalar) { if (scalar != 0) fixed (complex* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetSub(this complex[,] array, double scalar) { if (scalar != 0) fixed (complex* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetMul(this complex[,] array, double scalar) { if (scalar != 1) fixed (complex* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetDiv(this complex[,] array, double scalar) { if (scalar != 1) fixed (complex* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] LetConjMul(this complex[,] array, double scalar) { fixed (complex* p = array) Us.Cul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[,] Tocomplex(double[,] array0, double[,] array1) { Ex.SizeCheck(array0, array1); var result = array0.ZeroTo<double, complex>(); fixed (complex* r = result) fixed (double* p = array0, q = array1) Us.Tocomplex(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Tocomplex(this double[,] array) { var result = array.ZeroTo<double, complex>(); fixed (complex* r = result) fixed (double* p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static double[,] Real(this complex[,] array) { var result = array.ZeroTo<complex, double>(); fixed (double* r = result) fixed (complex* p = array) Us.Real(r, p, result.Length); return result; }
        public unsafe static double[,] Imag(this complex[,] array) { var result = array.ZeroTo<complex, double>(); fixed (double* r = result) fixed (complex* p = array) Us.Imag(r, p, result.Length); return result; }
        public unsafe static complex[,] Conj(this complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Conj(r, p, result.Length); return result; }
        public unsafe static complex[,] Pos(this complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static complex[,] Neg(this complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Neg(r, p, result.Length); return result; }
        public unsafe static complex[,] Rev(this complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Rev(r, p, result.Length); return result; }
        public unsafe static complex[,] Add(this complex[,] array0, complex[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Sub(this complex[,] array0, complex[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Mul(this complex[,] array0, complex[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Div(this complex[,] array0, complex[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static complex[,] ConjMul(this complex[,] array0, complex[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Cul(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Add(this complex[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Sub(this complex[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Mul(this complex[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Div(this complex[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static complex[,] ConjMul(this complex[,] array0, double[,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Cul(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Add(this double[,] array0, complex[,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Sub(this double[,] array0, complex[,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Sub(r, q, p, result.Length); return result; }
        public unsafe static complex[,] Mul(this double[,] array0, complex[,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[,] Div(this double[,] array0, complex[,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Div(r, q, p, result.Length); return result; }
        public unsafe static complex[,] Add(this complex[,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Sub(this complex[,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Mul(this complex[,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Div(this complex[,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] ConjMul(this complex[,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Cul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Add(complex scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Sub(complex scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static complex[,] Mul(complex scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Div(complex scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }
        public unsafe static complex[,] Add(this complex[,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Sub(this complex[,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Mul(this complex[,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Div(this complex[,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] ConjMul(this complex[,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Cul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Add(double scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Sub(double scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static complex[,] Mul(double scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[,] Div(double scalar, complex[,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }

        public unsafe static complex[,] LetAddMul(this complex[,] array0, complex[,] array1, double scalar) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.LetAddMul(p, q, scalar, array0.Length); return array0; }
        public unsafe static complex[,] LetMulAdd(this complex[,] array0, double scalar, complex[,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.LetMulAdd(p, scalar, q, array0.Length); return array0; }
        public unsafe static complex[,] AddMul(this complex[,] array0, complex[,] array1, double scalar) { Ex.SizeCheck(array0, array1); var result = array0.CloneX(); fixed (complex* r = result, p = array1) Us.LetAddMul(r, p, scalar, result.Length); return result; }
        #endregion

        #region complex[, ,]
        public unsafe static complex Sum(this complex[, ,] array) { fixed (complex* p = array) return Us.Sum(p, array.Length); }
        public unsafe static complex Average(this complex[, ,] array) { return array.Sum() / array.Length; }
        public unsafe static complex Inner(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static complex InnerConj(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumCul(p, q, array0.Length); }
        public unsafe static complex Inner(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) return Us.SumMul(p, q, array0.Length); }
        public unsafe static complex InnerConj(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) return Us.SumCul(p, q, array0.Length); }
        public unsafe static double Norm(this complex[, ,] array, double nu) { fixed (complex* p = array) return Math.Pow(Us.SumPowAbs(p, array.Length, nu), 1 / nu); }
        public unsafe static double Norm1(this complex[, ,] array) { fixed (complex* p = array) return Us.SumAbs(p, array.Length); }
        public unsafe static double Norm2(this complex[, ,] array) { return Math.Sqrt(SqNorm2(array)); }
        public unsafe static double SqNorm2(this complex[, ,] array) { fixed (complex* p = array) return Us.SumSqAbs(p, array.Length); }
        public unsafe static double MaxAbs(this complex[, ,] array) { return Math.Sqrt(MaxSqAbs(array)); }
        public unsafe static double MaxSqAbs(this complex[, ,] array) { fixed (complex* p = array) return Us.MaxSqAbs(p, array.Length); }
        public unsafe static double NormSub(complex[, ,] array0, complex[, ,] array1, double nu) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumPowSub(p, q, array0.Length, nu); }
        public unsafe static double Norm1Sub(complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumAbsSub(p, q, array0.Length); }
        public unsafe static double Norm2Sub(complex[, ,] array0, complex[, ,] array1) { return Math.Sqrt(SqNorm2Sub(array0, array1)); }
        public unsafe static double SqNorm2Sub(complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) return Us.SumSqSub(p, q, array0.Length); }

        public unsafe static complex[, ,] LetConj(this complex[, ,] array) { fixed (complex* p = array) Us.Conj(p, p, array.Length); return array; }
        public unsafe static complex[, ,] LetNeg(this complex[, ,] array) { fixed (complex* p = array) Us.Neg(p, p, array.Length); return array; }
        public unsafe static complex[, ,] LetRev(this complex[, ,] array) { fixed (complex* p = array) Us.LetRev(p, array.Length); return array; }
        public unsafe static complex[, ,] Let(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetAdd(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetSub(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetMul(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetDiv(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetSubr(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static complex[, ,] LetDivr(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static complex[, ,] LetConjMul(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Cul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetMulConj(this complex[, ,] array0, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.Cul(p, q, p, array0.Length); return array0; }
        public unsafe static complex[, ,] Let(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Let(p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetAdd(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Add(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetSub(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Sub(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetMul(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Mul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetDiv(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Div(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetConjMul(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Cul(p, p, q, array0.Length); return array0; }
        public unsafe static complex[, ,] LetSubr(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Sub(p, q, p, array0.Length); return array0; }
        public unsafe static complex[, ,] LetDivr(this complex[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0) fixed (double* q = array1) Us.Div(p, q, p, array0.Length); return array0; }
        public unsafe static complex[, ,] LetAdd(this complex[, ,] array, complex scalar) { if (scalar != 0) fixed (complex* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetSub(this complex[, ,] array, complex scalar) { if (scalar != 0) fixed (complex* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetMul(this complex[, ,] array, complex scalar) { if (scalar != 1) fixed (complex* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetDiv(this complex[, ,] array, complex scalar) { if (scalar != 1) fixed (complex* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetConjMul(this complex[, ,] array, complex scalar) { fixed (complex* p = array) Us.Cul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetSubr(this complex[, ,] array, complex scalar) { fixed (complex* p = array) Us.Sub(p, scalar, p, array.Length); return array; }
        public unsafe static complex[, ,] LetDivr(this complex[, ,] array, complex scalar) { fixed (complex* p = array) Us.Div(p, scalar, p, array.Length); return array; }
        public unsafe static complex[, ,] LetAdd(this complex[, ,] array, double scalar) { if (scalar != 0) fixed (complex* p = array) Us.Add(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetSub(this complex[, ,] array, double scalar) { if (scalar != 0) fixed (complex* p = array) Us.Sub(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetMul(this complex[, ,] array, double scalar) { if (scalar != 1) fixed (complex* p = array) Us.Mul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetDiv(this complex[, ,] array, double scalar) { if (scalar != 1) fixed (complex* p = array) Us.Div(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] LetConjMul(this complex[, ,] array, double scalar) { fixed (complex* p = array) Us.Cul(p, p, scalar, array.Length); return array; }
        public unsafe static complex[, ,] Tocomplex(double[, ,] array0, double[, ,] array1) { Ex.SizeCheck(array0, array1); var result = array0.ZeroTo<double, complex>(); fixed (complex* r = result) fixed (double* p = array0, q = array1) Us.Tocomplex(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Tocomplex(this double[, ,] array) { var result = array.ZeroTo<double, complex>(); fixed (complex* r = result) fixed (double* p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static double[, ,] Real(this complex[, ,] array) { var result = array.ZeroTo<complex, double>(); fixed (double* r = result) fixed (complex* p = array) Us.Real(r, p, result.Length); return result; }
        public unsafe static double[, ,] Imag(this complex[, ,] array) { var result = array.ZeroTo<complex, double>(); fixed (double* r = result) fixed (complex* p = array) Us.Imag(r, p, result.Length); return result; }
        public unsafe static complex[, ,] Conj(this complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Conj(r, p, result.Length); return result; }
        public unsafe static complex[, ,] Pos(this complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Let(r, p, result.Length); return result; }
        public unsafe static complex[, ,] Neg(this complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Neg(r, p, result.Length); return result; }
        public unsafe static complex[, ,] Rev(this complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Rev(r, p, result.Length); return result; }
        public unsafe static complex[, ,] Add(this complex[, ,] array0, complex[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Sub(this complex[, ,] array0, complex[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Mul(this complex[, ,] array0, complex[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Div(this complex[, ,] array0, complex[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] ConjMul(this complex[, ,] array0, complex[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0, q = array1) Us.Cul(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Add(this complex[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Sub(this complex[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Sub(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Mul(this complex[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Div(this complex[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Div(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] ConjMul(this complex[, ,] array0, double[, ,] array1) { var result = ToZero(array0, array1); fixed (complex* r = result, p = array0) fixed (double* q = array1) Us.Cul(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Add(this double[, ,] array0, complex[, ,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Add(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Sub(this double[, ,] array0, complex[, ,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Sub(r, q, p, result.Length); return result; }
        public unsafe static complex[, ,] Mul(this double[, ,] array0, complex[, ,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Mul(r, p, q, result.Length); return result; }
        public unsafe static complex[, ,] Div(this double[, ,] array0, complex[, ,] array1) { var result = ToZero(array1, array0); fixed (complex* r = result, p = array1) fixed (double* q = array0) Us.Div(r, q, p, result.Length); return result; }
        public unsafe static complex[, ,] Add(this complex[, ,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Sub(this complex[, ,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Mul(this complex[, ,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Div(this complex[, ,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] ConjMul(this complex[, ,] array, complex scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Cul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Add(complex scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Sub(complex scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static complex[, ,] Mul(complex scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Div(complex scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }
        public unsafe static complex[, ,] Add(this complex[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Sub(this complex[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Mul(this complex[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Div(this complex[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] ConjMul(this complex[, ,] array, double scalar) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Cul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Add(double scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Add(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Sub(double scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Sub(r, scalar, p, result.Length); return result; }
        public unsafe static complex[, ,] Mul(double scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Mul(r, p, scalar, result.Length); return result; }
        public unsafe static complex[, ,] Div(double scalar, complex[, ,] array) { var result = array.ZeroTo(); fixed (complex* r = result, p = array) Us.Div(r, scalar, p, result.Length); return result; }

        public unsafe static complex[, ,] LetAddMul(this complex[, ,] array0, complex[, ,] array1, double scalar) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.LetAddMul(p, q, scalar, array0.Length); return array0; }
        public unsafe static complex[, ,] LetMulAdd(this complex[, ,] array0, double scalar, complex[, ,] array1) { Ex.SizeCheck(array0, array1); fixed (complex* p = array0, q = array1) Us.LetMulAdd(p, scalar, q, array0.Length); return array0; }
        public unsafe static complex[, ,] AddMul(this complex[, ,] array0, complex[, ,] array1, double scalar) { Ex.SizeCheck(array0, array1); var result = array0.CloneX(); fixed (complex* r = result, p = array1) Us.LetAddMul(r, p, scalar, result.Length); return result; }
        #endregion

        #region double matrix
        public static double Inner(this double[,] matrix, double[] vector) { return Inner(vector, Multiply(matrix, vector)); }
        public static double Inner(this double[,] matrix, double[] vector0, double[] vector1) { return Inner(vector0, Multiply(matrix, vector1)); }

        public static double[,] AddDiag(this double[,] matrix, double[] vector) { return matrix.CloneX().LetAddDiag(vector); }
        public static double[,] LetAddDiag(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            if (matrix.GetLength(0) != n || matrix.GetLength(1) != n) ThrowException.SizeMismatch();
            for (int i = 0; i < n; i++) matrix[i, i] += vector[i];
            return matrix;
        }
        public static double[,] SubDiag(this double[,] matrix, double[] vector) { return matrix.CloneX().LetSubDiag(vector); }
        public static double[,] LetSubDiag(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            if (matrix.GetLength(0) != n || matrix.GetLength(1) != n) ThrowException.SizeMismatch();
            for (int i = 0; i < n; i++) matrix[i, i] -= vector[i];
            return matrix;
        }
        public static double[,] DiagSub(this double[] vector, double[,] matrix) { return matrix.CloneX().LetSubDiagLeft(vector); }
        public static double[,] SubDiagLeft(this double[,] matrix, double[] vector) { return matrix.CloneX().LetSubDiagLeft(vector); }
        public static double[,] LetSubDiagLeft(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            if (matrix.GetLength(0) != n || matrix.GetLength(1) != n) ThrowException.SizeMismatch();
            for (int i = 0; i < n; i++) matrix[i, i] = vector[i] - matrix[i, i];
            return matrix;
        }
        public static double[,] MulDiag(this double[,] matrix, double[] vector) { return matrix.CloneX().LetMulDiag(vector); }
        public static double[,] LetMulDiag(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            int m = matrix.GetLength(0);
            if (matrix.GetLength(1) != n) ThrowException.SizeMismatch();
            for (int i0 = 0; i0 < m; i0++)
                for (int i1 = 0; i1 < n; i1++)
                    matrix[i0, i1] *= vector[i1];
            return matrix;
        }
        public static double[,] DiagMul(this double[] vector, double[,] matrix) { return matrix.CloneX().LetMulDiagLeft(vector); }
        public static double[,] MulDiagLeft(this double[,] matrix, double[] vector) { return matrix.CloneX().LetMulDiagLeft(vector); }
        public static double[,] LetMulDiagLeft(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            int m = matrix.GetLength(1);
            if (matrix.GetLength(0) != n) ThrowException.SizeMismatch();
            for (int i0 = 0; i0 < n; i0++)
            {
                var a = vector[i0];
                for (int i1 = 0; i1 < m; i1++)
                    matrix[i0, i1] *= a;
            }
            return matrix;
        }
        public static double[,] DivDiag(this double[,] matrix, double[] vector) { return matrix.CloneX().LetDivDiag(vector); }
        public static double[,] LetDivDiag(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            int m = matrix.GetLength(0);
            if (matrix.GetLength(1) != n) ThrowException.SizeMismatch();
            for (int i0 = 0; i0 < m; i0++)
                for (int i1 = 0; i1 < n; i1++)
                    matrix[i0, i1] /= vector[i1];
            return matrix;
        }
        public static double[,] DiagDiv(this double[] vector, double[,] matrix) { return matrix.CloneX().LetDivDiagLeft(vector); }
        public static double[,] DivDiagLeft(this double[,] matrix, double[] vector) { return matrix.CloneX().LetDivDiagLeft(vector); }
        public static double[,] LetDivDiagLeft(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            int m = matrix.GetLength(1);
            if (matrix.GetLength(0) != n) ThrowException.SizeMismatch();
            for (int i0 = 0; i0 < n; i0++)
            {
                var a = 1 / vector[i0];
                for (int i1 = 0; i1 < m; i1++)
                    matrix[i0, i1] *= a;
            }
            return matrix;
        }

        public static double[,] Multiply(this double[] vector) { return Multiply(vector, vector); }
        public unsafe static double[,] Multiply(this double[] vector0, double[] vector1)
        {
            int n = vector1.Length;
            var result = new double[vector0.Length, n];
            fixed (double* r = result, p = vector0, q = vector1)
                for (int i = vector0.Length; --i >= 0; )
                    Us.Mul(&r[n * i], q, p[i], n);
            return result;
        }
        public unsafe static double[] Multiply(this double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            if (matrix.GetLength(1) != n) ThrowException.SizeMismatch();
            var result = new double[matrix.GetLength(0)];
            fixed (double* r = result, p = matrix, q = vector)
                for (int i = result.Length; --i >= 0; )
                    r[i] = Us.SumMul(&p[n * i], q, n);
            return result;
        }
        public unsafe static double[] Multiply(this double[] vector, double[,] matrix)
        {
            int n = matrix.GetLength(1);
            if (matrix.GetLength(0) != vector.Length) ThrowException.SizeMismatch();
            var result = new double[n];
            fixed (double* r = result, p = matrix, q = vector)
                for (int i = vector.Length; --i >= 0; )
                    Us.LetAddMul(r, &p[n * i], q[i], n);
            return result;
        }
        public unsafe static double[,] Multiply(this double[,] matrix0, double[,] matrix1)
        {
            int o = matrix0.GetLength(1);
            if (o != matrix1.GetLength(0)) ThrowException.SizeMismatch();
            int n = matrix0.GetLength(0);
            int m = matrix1.GetLength(1);
            var result = new double[n, m];
            fixed (double* r = result, p = matrix0, q = matrix1)
                for (int i = n; --i >= 0; )
                {
                    double* ri = &r[m * i];
                    double* pi = &p[o * i];
                    for (int j = o; --j >= 0; )
                        Us.LetAddMul(ri, &q[m * j], pi[j], m);
                }
            return result;
        }

        public static double Tr(this double[,] matrix)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1)) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix.GetLength(0); --i >= 0; ) a += matrix[i, i];
            return a;
        }
        public static double TrMultiply(double[,] matrix0, double[,] matrix1)
        {
            if (matrix0.GetLength(1) != matrix1.GetLength(0) || matrix0.GetLength(0) != matrix1.GetLength(1)) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix0.GetLength(0); --i >= 0; )
                for (int j = matrix0.GetLength(1); --j >= 0; )
                    a += matrix0[i, j] * matrix1[j, i];
            return a;
        }
        public static double TrMultiply(double[,] matrix0, double[,] matrix1, double[,] matrix2)
        {
            if (matrix0.GetLength(1) != matrix1.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(0) || matrix0.GetLength(0) != matrix2.GetLength(1)) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix1.GetLength(0); --i >= 0; )
                for (int j = matrix1.GetLength(1); --j >= 0; )
                {
                    double b = 0;
                    for (int k = matrix0.GetLength(0); --k >= 0; )
                        b += matrix0[k, i] * matrix2[j, k];
                    a += b * matrix1[i, j];
                }
            return a;
        }

        public static double[,] UnitMatrix(int size) { return DiagonalMatrix(size, 1); }
        public static double[,] DiagonalMatrix(int size, double value)
        {
            var result = new double[size, size];
            for (int i = size; --i >= 0; ) result[i, i] = value;
            return result;
        }
        public static double[,] DiagonalMatrix(double[] values)
        {
            var result = new double[values.Length, values.Length];
            for (int i = values.Length; --i >= 0; ) result[i, i] = values[i];
            return result;
        }
        public static void LetSymmetrical(this double[,] matrix)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1)) ThrowException.SizeMismatch();
            for (int i = matrix.GetLength(0); --i >= 0; )
                for (int j = i; --j >= 0; )
                {
                    double a = (matrix[i, j] + matrix[j, i]) * 0.5;
                    matrix[i, j] = a;
                    matrix[j, i] = a;
                }
        }
        public static double[,] T(this double[,] matrix)
        {
            var result = new double[matrix.GetLength(1), matrix.GetLength(0)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = matrix[j, i];
            return result;
        }

        // matrix inverse and determinant
        static double Det1by1(double[,] matrix)
        {
            return matrix[0, 0];
        }
        static double[,] Inverse1by1(double[,] matrix)
        {
            var det = matrix[0, 0];
            if (det == 0) ThrowException.WriteLine("Inverse: singular matrix");
            return new double[1, 1]{
                {1.0 / det}
            };
        }
        static double Det2by2(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1];
            double c = matrix[1, 0], d = matrix[1, 1];
            return a * d - b * c;
        }
        static double[,] Inverse2by2(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1];
            double c = matrix[1, 0], d = matrix[1, 1];
            var det = a * d - b * c;
            if (det == 0) ThrowException.WriteLine("Inverse: singular matrix");
            return new double[2, 2]{
                {d / det, -b / det},
                {-c / det, a / det}
            };
        }
        static double Det3by3(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1], c = matrix[0, 2];
            double d = matrix[1, 0], e = matrix[1, 1], f = matrix[1, 2];
            double g = matrix[2, 0], h = matrix[2, 1], i = matrix[2, 2];
            return a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g);
        }
        static double[,] Inverse3by3(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1], c = matrix[0, 2];
            double d = matrix[1, 0], e = matrix[1, 1], f = matrix[1, 2];
            double g = matrix[2, 0], h = matrix[2, 1], i = matrix[2, 2];
            double eifh = e * i - f * h, fgdi = f * g - d * i, dheg = d * h - e * g;
            var det = a * eifh + b * fgdi + c * dheg;
            if (det == 0) ThrowException.WriteLine("Inverse: singular matrix");
            return new double[3, 3]{
                { eifh / det, (c * h - b * i) / det, (b * f - c * e) / det},
                { fgdi / det, (a * i - c * g) / det, (c * d - a * f) / det},
                { dheg / det, (b * g - a * h) / det, (a * e - b * d) / det}
            };
        }
        static int[] LUDecomposition(double[][] matrix)
        {
            int n = matrix.GetLength(0);
            var pivot = new int[n];
            var vv = new double[n];
            for (int i = 0; i < n; i++)
            {
                //var max = Ex.Range(n).Max(j => Math.Abs(matrix[i][j]));
                var max = matrix[i].Max(v => Math.Abs(v));
                if (max == 0) ThrowException.WriteLine("LUDecomposition: singular matrix");
                vv[i] = 1 / max;
            }

            for (int k = 0; k < n; k++)
            {
                int p = k + Ex.FromTo(k, n).Select(i => vv[i] * Math.Abs(matrix[i][k])).MaxIndex();
                if (p != k)
                {
                    Mt.Swap(ref matrix[p], ref matrix[k]);
                    vv[p] = vv[k];
                }
                pivot[k] = p;
                var mk = matrix[k];
                if (mk[k] == 0) { ThrowException.WriteLine("LUDecomposition: singular matrix"); mk[k] = 1e-40; }
                for (int i = k; ++i < n; )  //順不同並列可
                {
                    var mi = matrix[i];
                    var temp = mi[k] /= mk[k];
                    for (int j = k; ++j < n; )
                        mi[j] -= temp * mk[j];
                }
            }
            return pivot;
        }
        public static Tuple<double[][], int[]> LUDecomposition(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n != matrix.GetLength(1)) ThrowException.SizeMismatch();
            double[][] LU = New.Array(n, i => New.Array(n, j => matrix[i, j]));
            int[] pivot = LUDecomposition(LU);
            return Tuple.Create(LU, pivot);
        }
        public static double[] Divide(this double[] vector, Tuple<double[][], int[]> LUinfo)
        {
            double[][] LU = LUinfo.Item1;
            int[] pivot = LUinfo.Item2;
            int n = vector.Length;
            if (n != pivot.Length) ThrowException.SizeMismatch();

            var r = vector.CloneX();
            for (int ii = 0, i = 0; i < n; i++)
            {
                var sum = r[pivot[i]];
                r[pivot[i]] = r[i];
                if (ii != 0)
                    for (int j = ii - 1; j < i; j++) sum -= LU[i][j] * r[j];
                else if (sum != 0)
                    ii = i + 1;
                r[i] = sum;
            }
            for (int i = n - 1; i >= 0; i--)
            {
                var sum = r[i];
                for (int j = i + 1; j < n; j++) sum -= LU[i][j] * r[j];
                r[i] = sum / LU[i][i];
            }
            return r;
        }
        public static double[] Divide(this double[] vector, double[,] matrix) { return vector.Divide(LUDecomposition(matrix)); }

        public static double[,] Inverse(this double[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n != matrix.GetLength(1)) ThrowException.SizeMismatch();
            if (n == 1) return Inverse1by1(matrix);
            if (n == 2) return Inverse2by2(matrix);
            if (n == 3) return Inverse3by3(matrix);

            double[][] LU = New.Array(n, i => New.Array(n, j => matrix[j, i]));
            int[] pivot = LUDecomposition(LU);
            int[] index = New.Array(n, i => i);
            for (int i = 0; i < n; i++)
                Mt.Swap(ref index[i], ref index[pivot[i]]);

            var result = new double[n, n];
            var vec = new double[n];
            for (int i = n; --i >= 0; )  //順不同並列可
            {
                vec.Clear();
                vec[i] = 1.0;
                for (int j = i; ++j < n; )
                {
                    double sum = 0;
                    for (int k = i; k < j; k++)
                        sum -= LU[j][k] * vec[k];
                    vec[j] = sum;
                }
                for (int j = n; --j >= 0; )
                {
                    double sum = vec[j];
                    for (int k = j; ++k < n; )
                        sum -= LU[j][k] * vec[k];
                    vec[j] = sum / LU[j][j];
                }
                int idx = index[i];
                for (int j = n; --j >= 0; )
                    result[idx, j] = vec[j];
            }
            return result;
        }

        public static double Det(this double[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n == 1) return Det1by1(matrix);
            if (n == 2) return Det2by2(matrix);
            if (n == 3) return Det3by3(matrix);

            double[][] L = New.Array(n, i => New.Array(n, j => matrix[i, j]));
            int[] pivot = LUDecomposition(L);
            int swap = Ex.Range(n).Count(i => pivot[i] != i);
            return (swap % 2 == 0 ? 1 : -1) * Mt.Product(n, i => L[i][i]);
        }
        public static double LogDet(this double[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n <= 3) return Math.Log(Mt.Det(matrix));

            double[][] L = New.Array(n, i => New.Array(n, j => matrix[i, j]));
            int[] pivot = LUDecomposition(L);
            int swap = Ex.Range(n).Count(i => pivot[i] != i);
            double sign = swap % 2 == 0 ? 1.0 : -1.0;
            if (sign == -1) ThrowException.Argument("LogDet: determinant is not positive");
            return Mt.Sum(n, i => Math.Log(L[i][i]));
        }

        public static double[,] PseudoInverse(this double[,] matrix)
        {
            var transpose = matrix.T();
            if (matrix.GetLength(0) >= matrix.GetLength(1))
                return Mt.Multiply(Mt.Multiply(transpose, matrix).Inverse(), transpose);
            else
                return Mt.Multiply(transpose, Mt.Multiply(matrix, transpose).Inverse());
        }

        // Inverse of symmetric positive definite matrix
        static double[,] CholeskyDecomposition(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            var L = new double[n, n];
            //対角と下三角部分を計算
            for (int j = 0; j < n; j++)
            {
                for (int i = j; i < n; i++)
                {
                    double sum = matrix[j, i];
                    for (int k = j; --k >= 0; )
                        sum -= L[j, k] * L[i, k];
                    if (j == i && sum <= 0.0) ThrowException.WriteLine("CholeskyDecomposition: not a positive definite matrix");
                    L[i, j] = (j == i) ? Math.Sqrt(sum) : sum / L[j, j];
                }
            }
            return L;
        }
        public static double[,] InverseSymmetricPositiveDefinite(this double[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n != matrix.GetLength(1)) ThrowException.SizeMismatch();
            double[,] L = CholeskyDecomposition(matrix);
            //対角と上三角部分を計算
            for (int j = 0; j < n; j++)
            {
                L[j, j] = 1 / L[j, j];
                for (int i = j; --i >= 0; )
                {
                    double sum = 0;
                    for (int k = j; --k >= i; )
                        sum -= L[j, k] * L[i, k];
                    L[i, j] = sum * L[j, j];
                }
            }
            for (int j = n; --j >= 0; )
            {
                for (int i = 0; i <= j; i++)
                {
                    double sum = L[i, j];
                    for (int k = j; ++k < n; )
                        sum -= L[k, j] * L[i, k];
                    L[i, j] = sum * L[j, j];
                }
            }
            //下三角部分を上三角部分からコピー
            for (int i = n; --i >= 0; )
                for (int j = i; --j >= 0; )
                    L[i, j] = L[j, i];
            return L;
        }

        static void _QR(double[] vector0, double[] vector1, double c, double s)
        {
            for (int i = vector0.Length; --i >= 0; )
            {
                var x = vector0[i];
                var y = vector1[i];
                vector0[i] = x * c + y * s;
                vector1[i] = y * c - x * s;
            }
        }

        //M: symmetric matrix
        //Mを転置させた計算の方が速いため改変
        static double[] Householder(double[][] M, double[] D)
        {
            int n = M.Length;
            var E = new double[n];
            for (int i = n; --i >= 2; )
            {
                var Mi = new double[i];
                for (int l = i; --l >= 0; ) Mi[l] = M[l][i];
                double scale = Mi.Sum(x => Math.Abs(x));
                if (scale == 0) continue;
                for (int l = i; --l >= 0; ) Mi[l] /= scale;

                double Di;
                {
                    var f = Mi.Sum(x => x * x);
                    var g = Math.Sqrt(f);
                    if (Mi[i - 1] < 0) g *= -1;
                    E[i] = -g * scale;
                    D[i] = Di = f + g * Mi[i - 1];
                    Mi[i - 1] += g;
                }

                var hh = 0.0;
                for (int j = 0; j < i; j++)
                {
                    M[i][j] = Mi[j] / Di;
                    {
                        var a = 0.0;
                        for (int l = i; --l > j; ) a += M[j][l] * Mi[l];
                        for (int l = j + 1; --l >= 0; ) a += M[l][j] * Mi[l];
                        E[j] = a / Di;
                    }
                    hh += E[j] * Mi[j];
                }
                hh /= 2 * Di;
                for (int j = 0; j < i; j++)
                {
                    E[j] -= Mi[j] * hh;
                    for (int k = j + 1; --k >= 0; ) M[k][j] -= Mi[j] * E[k] + Mi[k] * E[j];
                }
                for (int l = i; --l >= 0; ) M[l][i] = Mi[l];
            }
            E[1] = M[0][1];

            for (int i = 0; i < n; i++)
            {
                if (D[i] != 0)
                {
                    for (int j = 0; j < i; j++)
                    {
                        var a = 0.0;
                        for (int l = i; --l >= 0; ) a += M[j][l] * M[l][i];
                        for (int l = i; --l >= 0; ) M[j][l] -= M[i][l] * a;
                    }
                }
                D[i] = M[i][i];
                M[i][i] = 1;
                for (int l = i; --l >= 0; ) M[i][l] = 0;
                for (int l = i; --l >= 0; ) M[l][i] = 0;
            }
            return E;
        }
        static void QLMethod(double[][] M, double[] D, double[] E)
        {
            int N = M.Length;
            for (int l = 0; l < N; l++)
            {
                while (true)
                {
                    int m;
                    for (m = l; m < N - 1; m++)
                        if (IsTooSmall(E[m], Math.Abs(D[m]) + Math.Abs(D[m + 1])) || Math.Abs(E[m]) < Mt.DoubleEpsilon) break;
                    if (m == l) break;

                    var g = (D[l + 1] - D[l]) / E[l] * 0.5;
                    g += g * Math.Sqrt(1 + 1 / (g * g));
                    g = E[l] / g + D[m] - D[l];
                    double s = 1, c = 1, p = 0;
                    int i;
                    for (i = m; --i >= l; )
                    {
                        double b = E[i] * c;
                        {
                            var f = E[i] * s;
                            var py = Norm2(f, g);
                            E[i + 1] = py;
                            if (py == 0) break;
                            s = f / py;
                            c = g / py;
                        }
                        {
                            var q = D[i + 1] - p;
                            var r = (D[i] - q) * s + b * c * 2;
                            p = r * s;
                            D[i + 1] = q + p;
                            g = r * c - b;
                        }
                        _QR(M[i], M[i + 1], c, -s);
                    }
                    if (i < l) E[l] = g;
                    D[i + 1] -= p;
                    E[m] = 0;
                }
            }
        }
        //matrix: symmetric matrix
        //matrix => V * diag[D] * V^T
        //固有値の配列と固有ベクトルの配列を返す
        public static Tuple<double[], double[][]> EigenValueDecomposition(double[,] matrix)
        {
            int N = matrix.GetLength(0);
            if (N != matrix.GetLength(1)) ThrowException.SizeMismatch();
            var D = new double[N];
            var V = New.Array(N, i => New.Array(N, j => matrix[i, j]));

            double[] E = Householder(V, D);
            for (int i = 0; i < N - 1; i++) E[i] = E[i + 1];
            E[N - 1] = 0;
            QLMethod(V, D, E);

            int[] order = D.SortIndex().LetReverse();
            D.LetSortAsIndex(order);
            V.LetSortAsIndex(order);
            for (int i = N; --i >= 0; )
                if (V[i].Count(a => a < 0) > N / 2) V[i].LetNeg();
            return Tuple.Create(D, V);
        }
        //matrix => V * diag[D] * V^T
        //V, Dを返す
        public static Tuple<double[,], double[]> EigenValueDecompositionM(double[,] matrix)
        {
            var T = EigenValueDecomposition(matrix);
            return Tuple.Create(
                New.Array(T.Item2[0].Length, T.Item2.Length, (i, j) => T.Item2[j][i]),
                T.Item1
            );
        }

        //matrix -> U diag[W] V^T
        //U, V: 縦ベクトルの配列を返す
        public static Tuple<double[][], double[], double[][]> SingularValueDecomposition(double[,] matrix)
        {
            int XN = matrix.GetLength(1);
            int YN = matrix.GetLength(0);
            int ZN = Math.Min(YN, XN);
            var U = New.Array(XN, i => New.Array(YN, j => matrix[j, i]));
            var V = New.Array(XN, i => new double[XN]);
            var W = new double[XN];
            var R = new double[XN];
            double anorm = 0;
            for (int i = 0; i < XN; i++)
            {
                var I = Ex.FromTo(i, YN);
                if (i < YN)
                {
                    double scale = I.Sum(k => Math.Abs(U[i][k]));
                    if (scale != 0)
                    {
                        foreach (var k in I) U[i][k] /= scale;
                        var s = I.Sum(k => Mt.Sq(U[i][k]));
                        var f = U[i][i];
                        var g = Math.Sqrt(s);
                        if (f > 0) g *= -1;
                        var h = f * g - s;
                        U[i][i] = f - g;
                        for (int j = i; ++j < XN; )
                        {
                            var coef = 0.0;
                            for (int k = i; k < YN; k++) coef += U[j][k] * U[i][k];
                            coef /= h;
                            for (int k = i; k < YN; k++) U[j][k] += coef * U[i][k];
                        }
                        foreach (var k in I) U[i][k] *= scale;
                        W[i] = scale * g;
                    }
                }
                int l = i + 1;
                var L = Ex.FromTo(l, XN);
                if (i < YN && i < XN - 1)
                {
                    double scale = L.Sum(k => Math.Abs(U[k][i]));
                    if (scale != 0.0)
                    {
                        foreach (var k in L) U[k][i] /= scale;
                        var s = L.Sum(k => Mt.Sq(U[k][i]));
                        var f = U[l][i];
                        var g = -Math.Sign(f) * Math.Sqrt(s);
                        var h = f * g - s;
                        U[l][i] = f - g;
                        foreach (var k in L) R[k] = U[k][i] / h;
                        for (int j = l; j < YN; j++)
                        {
                            double coef = 0;
                            for (int k = l; k < XN; k++) coef += U[k][j] * U[k][i];
                            for (int k = l; k < XN; k++) U[k][j] += coef * R[k];
                        }
                        foreach (var k in L) U[k][i] *= scale;
                        R[l] = scale * g;
                    }
                }
                Mt.LetMax(ref anorm, Math.Abs(W[i]) + Math.Abs(R[i]));
                if (i > YN) R[i] = 0.0;
            }

            for (int i = XN; --i >= 0; )
            {
                if (i < XN - 1)
                {
                    int l = i + 1;
                    var L = Ex.FromTo(l, XN);
                    var g = R[l];
                    if (g != 0.0)
                    {
                        if (U[l][i] != 0.0)
                        {
                            var coef = 1 / U[l][i] / g;
                            foreach (var k in L) V[i][k] = U[k][i] * coef;
                        }
                        for (int j = l; j < XN; j++)
                        {
                            var coef = 0.0;
                            for (int k = l; k < XN; k++) coef += V[j][k] * U[k][i];
                            for (int k = l; k < XN; k++) V[j][k] += coef * V[i][k];
                        }
                    }
                    foreach (var k in L) V[i][k] = 0;
                    foreach (var k in L) V[k][i] = 0;
                }
                V[i][i] = 1.0;
            }

            for (int i = ZN; --i >= 0; )
            {
                int l = i + 1;
                var I = Ex.FromTo(i, YN);
                foreach (var k in Ex.FromTo(l, XN)) U[k][i] = 0;
                var g = W[i];
                if (g != 0.0)
                {
                    g = 1 / g;
                    for (int j = l; j < XN; j++)
                    {
                        var coef = 0.0;
                        for (int k = l; k < YN; k++) coef += U[j][k] * U[i][k];
                        coef = coef / U[i][i] * g;
                        for (int k = i; k < YN; k++) U[j][k] += coef * U[i][k];  //
                    }
                    foreach (var k in I) U[i][k] *= g;
                }
                else
                    foreach (var k in I) U[i][k] = 0;
                U[i][i] += 1;
            }

            for (int k = XN; --k >= 0; )
            {
                for (int iterations = 30; ; iterations--)
                {
                    int l;
                    for (l = k; l > 0; l--)
                    {  // R[0]==0;
                        if (IsTooSmall(R[l], anorm)) break;
                        if (IsTooSmall(W[l - 1], anorm))
                        {
                            double c = 0, s = 1;
                            for (int i = l; i <= k; i++)
                            {
                                var f = s * R[i];
                                R[i] *= c;
                                if (IsTooSmall(f, anorm)) break;
                                var g = W[i];
                                var py = Norm2(f, g);
                                W[i] = py;
                                c = g / py;
                                s = -f / py;
                                _QR(U[l - 1], U[i], c, s);
                            }
                            break;
                        }
                    }

                    if (l == k)
                    {
                        if (W[k] < 0) { W[k] *= -1; Mt.LetMul(V[k], -1); }
                        break;
                    }
                    if (iterations == 0) { ThrowException.WriteLine("SingularValueDecomposition: too many iterations"); break; }

                    {
                        var f = 0.0;
                        var x = W[l];
                        {
                            var y = W[k - 1];
                            var z = W[k];
                            var g = R[k - 1];
                            var h = R[k];
                            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                            var py = f * (1 + Math.Sqrt(1 + 1 / (f * f)));
                            f = ((x - z) * (x + z) + h * (y / py - h)) / x;
                        }
                        double c = 1, s = 1;
                        for (int j = l; j < k; j++)
                        {
                            int i = j + 1;
                            var g = R[i] * c;
                            var h = R[i] * s;
                            {
                                var py = Norm2(f, h);
                                R[j] = py;
                                c = f / py;
                                s = h / py;
                            }
                            f = x * c + g * s;
                            g = g * c - x * s;

                            var y = W[i] * c;
                            h = W[i] * s;
                            _QR(V[j], V[i], c, s);
                            {
                                var py = Norm2(f, h);
                                W[j] = py;
                                if (py != 0)
                                {
                                    c = f / py;
                                    s = h / py;
                                }
                            }
                            f = g * c + y * s;
                            x = y * c - g * s;
                            _QR(U[j], U[i], c, s);
                        }
                        R[l] = 0;
                        R[k] = f;
                        W[k] = x;
                    }
                }
            }
            {
                var order = W.SortIndex().LetReverse();
                W.LetSortAsIndex(order);
                U.LetSortAsIndex(order);
                V.LetSortAsIndex(order);

                for (int k = 0; k < XN; k++)
                {
                    int s = U[k].Count(a => a < 0) + V[k].Count(a => a < 0);
                    if (s > (YN + XN) / 2)
                    {
                        Mt.LetMul(U[k], -1);
                        Mt.LetMul(V[k], -1);
                    }
                }
            }
            return Tuple.Create(U, W, V);
        }
        //matrix -> U diag[W] V^T
        //U, W, Vを返す
        public static Tuple<double[,], double[], double[,]> SingularValueDecompositionM(double[,] matrix)
        {
            var T = SingularValueDecomposition(matrix);
            return Tuple.Create(
                New.Array(T.Item1[0].Length, T.Item1.Length, (i, j) => T.Item1[j][i]),
                T.Item2,
                New.Array(T.Item3.Length, T.Item3.Length, (i, j) => T.Item3[j][i])
            );
        }

        // Rotation Matrix
        public static Double2 Rotate(Double2 vector, double angle)
        {
            var cos = Math.Cos(angle);
            var sin = Math.Sin(angle);
            return new Double2(cos * vector.X - sin * vector.Y, sin * vector.X + cos * vector.Y);
        }
        // 右手座標系ならvectorをaxisを軸にして右ねじの回転方向に回転させる．axisの長さ=1
        public static Double3 Rotate(Double3 vector, Double3 axis, double angle)
        {
            var cos = Math.Cos(angle);
            var sin = Math.Sin(angle);
            var cos1 = 1 - cos;
            var x = axis.X;
            var y = axis.Y;
            var z = axis.Z;
            return new Double3(
                (cos1 * x * x + cos) * vector.X + (cos1 * x * y - sin * z) * vector.Y + (cos1 * z * x + sin * y) * vector.Z,
                (cos1 * x * y + sin * z) * vector.X + (cos1 * y * y + cos) * vector.Y + (cos1 * y * z - sin * x) * vector.Z,
                (cos1 * z * x - sin * y) * vector.X + (cos1 * y * z + sin * x) * vector.Y + (cos1 * z * z + cos) * vector.Z
            );
        }
        public static double[,] RotationMatrix(double angle)
        {
            var cos = Math.Cos(angle);
            var sin = Math.Sin(angle);
            var matrix = new double[2, 2]{
                { cos, -sin },
                { sin, cos }
            };
            return matrix;
        }
        public static double[,] RotationMatrix(Double3 axis, double angle)
        {
            var cos = Math.Cos(angle);
            var sin = Math.Sin(angle);
            var cos1 = 1 - cos;
            var x = axis.X;
            var y = axis.Y;
            var z = axis.Z;
            var matrix = new double[3, 3] {
                { cos1 * x * x + cos, cos1 * x * y - sin * z, cos1 * z * x + sin * y },
                { cos1 * x * y + sin * z, cos1 * y * y + cos, cos1 * y * z - sin * x },
                { cos1 * z * x - sin * y, cos1 * y * z + sin * x, cos1 * z * z + cos },
            };
            return matrix;
        }
        #endregion

        #region complex matrix
        public static double[,] SqAbs(this complex[,] array)
        {
            var result = new double[array.GetLength(0), array.GetLength(1)];
            unsafe { fixed (double* r = result) fixed (complex* p = array) Us.SqAbs(r, p, result.Length); } return result;
        }

        public static complex[,] UnitMatrixcomplex(int size) { return DiagonalMatrix(size, complex.One); }
        public static complex[,] DiagonalMatrix(int size, complex value)
        {
            var result = new complex[size, size];
            for (int i = size; --i >= 0; ) result[i, i] = value;
            return result;
        }
        public static complex[,] DiagonalMatrix(complex[] values)
        {
            var result = new complex[values.Length, values.Length];
            for (int i = values.Length; --i >= 0; ) result[i, i] = values[i];
            return result;
        }
        public static void LetSymmetrical(this complex[,] matrix)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1)) ThrowException.SizeMismatch();
            for (int i = matrix.GetLength(0); --i >= 0; )
                for (int j = i; --j >= 0; )
                {
                    var a = (matrix[i, j].Real + matrix[j, i].Real) * 0.5;
                    var b = (matrix[i, j].Imaginary - matrix[j, i].Imaginary) * 0.5;
                    matrix[i, j] = new complex(a, +b);
                    matrix[j, i] = new complex(a, -b);
                }
        }
        public static complex[,] T(this complex[,] matrix)
        {
            var result = new complex[matrix.GetLength(1), matrix.GetLength(0)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = matrix[j, i];
            return result;
        }
        public static complex[,] H(this complex[,] matrix)
        {
            var result = new complex[matrix.GetLength(1), matrix.GetLength(0)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = complex.Conjugate(matrix[j, i]);
            return result;
        }

        public unsafe static complex[] Multiply(complex[,] matrix, complex[] vector)
        {
            int n = vector.Length;
            if (matrix.GetLength(1) != n) ThrowException.SizeMismatch();
            var result = new complex[matrix.GetLength(0)];
            fixed (complex* r = result, p = matrix, q = vector)
                for (int i = result.Length; --i >= 0; )
                    r[i] = Us.SumMul(&p[n * i], q, n);
            return result;
        }
        public unsafe static complex[] Multiply(complex[] vector, complex[,] matrix)
        {
            int n = matrix.GetLength(1);
            if (matrix.GetLength(0) != vector.Length) ThrowException.SizeMismatch();
            var result = new complex[n];
            fixed (complex* r = result, p = matrix, q = vector)
                for (int i = vector.Length; --i >= 0; )
                    Us.LetAddMul(r, &p[n * i], q[i], n);
            return result;
        }
        public unsafe static complex[,] Multiply(complex[,] matrix0, complex[,] matrix1)
        {
            int o = matrix0.GetLength(0);
            int n = matrix1.GetLength(0);
            int m = matrix1.GetLength(1);
            if (n != matrix0.GetLength(1)) ThrowException.SizeMismatch();
            var result = new complex[o, m];
            fixed (complex* r = result, p = matrix0, q = matrix1)
                for (int i = o; --i >= 0; )
                {
                    complex* ri = &r[m * i];
                    complex* pi = &p[n * i];
                    for (int j = n; --j >= 0; )
                        Us.LetAddMul(ri, &q[m * j], pi[j], m);
                }
            return result;
        }

        static int[] LUDecomposition(complex[][] matrix)
        {
            var n = matrix.GetLength(0);
            var pivot = new int[n];
            var vv = new complex[n];
            for (int i = 0; i < n; i++)
            {
                double max = Ex.Range(n).Max(j => complex.Abs(matrix[i][j]));
                if (max == 0) ThrowException.WriteLine("LUDecomposition: singular matrix");
                vv[i] = (complex)(1 / max);
            }

            for (int k = 0; k < n; k++)
            {
                int p = k + Ex.FromTo(k, n).Select(i => complex.Abs(vv[i] * matrix[i][k])).MaxIndex();
                if (p != k)
                {
                    Mt.Swap(ref matrix[p], ref matrix[k]);
                    vv[p] = vv[k];
                }
                pivot[k] = p;
                var mk = matrix[k];
                if (mk[k] == 0.0) { ThrowException.WriteLine("LUDecomposition: singular matrix"); mk[k] = (complex)1e-40; }
                for (int i = k; ++i < n; )  //順不同並列可
                {
                    var mi = matrix[i];
                    var temp = mi[k] /= mk[k];
                    for (int j = k; ++j < n; )
                        mi[j] -= temp * mk[j];
                }
            }
            return pivot;
        }
        public static complex[,] Inverse(this complex[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n != matrix.GetLength(1)) ThrowException.SizeMismatch();
            //if (n == 1) return Inverse1by1(matrix);
            //if (n == 2) return Inverse2by2(matrix);
            //if (n == 3) return Inverse3by3(matrix);

            complex[][] LU = New.Array(n, i => New.Array(n, j => matrix[j, i]));
            int[] pivot = LUDecomposition(LU);
            int[] index = New.Array(n, i => i);
            for (int i = 0; i < n; i++)
                Mt.Swap(ref index[i], ref index[pivot[i]]);

            var result = new complex[n, n];
            var vec = new complex[n];
            for (int i = n; --i >= 0; )  //順不同並列可
            {
                vec.Clear();
                vec[i] = complex.One;
                for (int j = i; ++j < n; )
                {
                    var sum = complex.Zero;
                    for (int k = i; k < j; k++)
                        sum -= LU[j][k] * vec[k];
                    vec[j] = sum;
                }
                for (int j = n; --j >= 0; )
                {
                    var sum = vec[j];
                    for (int k = j; ++k < n; )
                        sum -= LU[j][k] * vec[k];
                    vec[j] = sum / LU[j][j];
                }
                int idx = index[i];
                for (int j = n; --j >= 0; )
                    result[idx, j] = vec[j];
            }
            return result;
        }

        public static complex[,] PseudoInverse(this complex[,] matrix)
        {
            var transpose = matrix.H();
            if (matrix.GetLength(0) >= matrix.GetLength(1))
                return Mt.Multiply(Mt.Multiply(transpose, matrix).Inverse(), transpose);
            else
                return Mt.Multiply(transpose, Mt.Multiply(matrix, transpose).Inverse());
        }
        #endregion



        #endregion
    }

    // Numerical functions
    public static partial class Nm
    {
        #region Spline classes
        public class Spline3
        {
            double[] DataX, DataY, DataZ;
            public Spline3() { }
            public Spline3(IEnumerable<double> dataX, IEnumerable<double> dataY) { Set(dataX, dataY); }

            int FindSection(double value)
            {
                int i = 0;
                for (int j = DataX.Length - 1; i < j; )
                {
                    int k = (i + j) / 2;
                    if (DataX[k] < value) i = k + 1; else j = k;
                }
                if (i > 0) --i;
                return i;
            }
            public double Interpolate(double value)
            {
                int i = FindSection(value);
                double h = DataX[i + 1] - DataX[i];
                double d = value - DataX[i];
                return (((DataZ[i + 1] - DataZ[i]) * d / h + DataZ[i] * 3) * d + ((DataY[i + 1] - DataY[i]) / h - (DataZ[i] * 2 + DataZ[i + 1]) * h)) * d + DataY[i];
            }
            public double CalcGrad(double value)
            {
                int i = FindSection(value);
                double h = DataX[i + 1] - DataX[i];
                double d = value - DataX[i];
                return ((DataZ[i + 1] - DataZ[i]) * d / h * 3 + DataZ[i] * 6) * d + ((DataY[i + 1] - DataY[i]) / h - (DataZ[i] * 2 + DataZ[i + 1]) * h);
            }

            public void Set(IEnumerable<double> dataX, IEnumerable<double> dataY)
            {
                if (dataY == null) ThrowException.Argument("dataY");
                DataY = dataY.ToArray();
                int N = DataY.Length;
                if (N == 0) ThrowException.Argument("dataY");
                if (N == 1) DataY = new double[] { DataY[0], DataY[0], DataY[0] };
                if (N == 2) DataY = new double[] { DataY[0], (DataY[0] + DataY[1]) / 2, DataY[1] };

                if (dataX != null)
                {
                    DataX = dataX.ToArray();
                    if (DataX.Length == 1) DataX = new double[] { DataX[0], DataX[0], DataX[0] };
                    if (DataX.Length == 2) DataX = new double[] { DataX[0], (DataX[0] + DataX[1]) / 2, DataX[1] };
                }
                if (dataX == null || DataX.Length == 0)
                {
                    DataX = new double[N];
                    DataX[0] = 0;
                    for (int i = 1; i < N; i++) DataX[i] = DataX[i - 1] + Math.Max(Math.Abs(DataY[i] - DataY[i - 1]), 1e-10);
                    for (int i = 1; i < N; i++) DataX[i] /= DataX[N - 1];
                }
                if (DataX.Length != N) ThrowException.ArgumentOutOfRange("dataX, dataY");

                DataZ = new double[N];
                double[] h = new double[N];
                double[] d = new double[N];
                DataZ[0] = DataZ[N - 1] = 0;
                for (int i = 0; i < N - 1; i++)
                {
                    h[i] = DataX[i + 1] - DataX[i];
                    d[i + 1] = (DataY[i + 1] - DataY[i]) / h[i];
                }
                DataZ[1] = d[2] - d[1] - h[0] * DataZ[0];
                d[1] = 2 * (DataX[2] - DataX[0]);
                for (int i = 1; i < N - 2; i++)
                {
                    double t = h[i] / d[i];
                    DataZ[i + 1] = d[i + 2] - d[i + 1] - DataZ[i] * t;
                    d[i + 1] = 2 * (DataX[i + 2] - DataX[i]) - h[i] * t;
                }
                DataZ[N - 2] -= h[N - 2] * DataZ[N - 1];
                for (int i = N - 2; i > 0; i--)
                    DataZ[i] = (DataZ[i] - h[i] * DataZ[i + 1]) / d[i];
            }
        }

        public class Spline1
        {
            double[] DataX, DataY;
            public Spline1() { }
            public Spline1(IEnumerable<double> dataX, IEnumerable<double> dataY) { Set(dataX, dataY); }

            int FindSection(double value)
            {
                int i = 0;
                for (int j = DataX.Length - 1; i < j; )
                {
                    int k = (i + j) / 2;
                    if (DataX[k] < value) i = k + 1; else j = k;
                }
                if (i > 0) --i;
                return i;
            }
            public double Interpolate(double value)
            {
                int i = FindSection(value);
                return (DataY[i + 1] - DataY[i]) / (DataX[i + 1] - DataX[i]) * (value - DataX[i]) + DataY[i];
            }
            public void Set(IEnumerable<double> dataX, IEnumerable<double> dataY)
            {
                if (dataX == null || dataY == null) ThrowException.Argument("dataX, dataY");
                DataX = dataX.ToArray();
                DataY = dataY.ToArray();
                if (DataX.Length == 0 || DataY.Length == 0 || DataX.Length != DataY.Length) ThrowException.Argument("dataX, dataY");
            }
        }
        #endregion

        #region Interporation
        // NumericalRecipes3.Base_Interp
        public abstract class InterporationBase
        {
            protected int n;
            protected int mm;
            protected int jsav = 0;
            protected int cor = 0;
            protected int dj;
            protected double[] xx;
            protected double[] yy;

            protected InterporationBase(double[] x, double[] y, int m)
            {
                this.n = x.Length;
                this.mm = m;
                this.xx = x;
                this.yy = y;
                this.dj = Math.Min(1, (int)Math.Pow((double)n, 0.25));
            }
            protected double interp(double x)
            {
                int jlo = cor != 0 ? hunt(x) : locate(x);
                return rawinterp(jlo, x);
            }
            public int locate(double x)
            {
                if (n < 2 || mm < 2 || mm > n) { ThrowException.Argument("InterporationBase: locate size error"); return 0; }
                bool ascnd = (xx[n - 1] >= xx[0]);
                var jl = 0;
                var ju = n - 1;
                while (ju - jl > 1)
                {
                    var jm = (ju + jl) >> 1;
                    if (x >= xx[jm] == ascnd) jl = jm;
                    else ju = jm;
                }
                cor = Math.Abs(jl - jsav) > dj ? 0 : 1;
                jsav = jl;
                return Math.Max(0, Math.Min(n - mm, jl - ((mm - 2) >> 1)));
            }
            public int hunt(double x)
            {
                int jl = jsav, jm, ju, inc = 1;
                if (n < 2 || mm < 2 || mm > n) { ThrowException.Argument("InterporationBase: hunt size error"); return 0; }
                bool ascnd = (xx[n - 1] >= xx[0]);
                if (jl < 0 || jl > n - 1) { jl = 0; ju = n - 1; }
                else
                {
                    if (x >= xx[jl] == ascnd)
                    {
                        for (; ; )
                        {
                            ju = jl + inc;
                            if (ju >= n - 1) { ju = n - 1; break; }
                            if (x < xx[ju] == ascnd) break;
                            jl = ju; inc += inc;
                        }
                    }
                    else
                    {
                        ju = jl;
                        for (; ; )
                        {
                            jl = jl - inc;
                            if (jl <= 0) { jl = 0; break; }
                            if (x >= xx[jl] == ascnd) break;
                            ju = jl; inc += inc;
                        }
                    }
                }
                while (ju - jl > 1)
                {
                    jm = (ju + jl) >> 1;
                    if (x >= xx[jm] == ascnd) jl = jm;
                    else ju = jm;
                }
                cor = Math.Abs(jl - jsav) > dj ? 0 : 1;
                jsav = jl;
                return Math.Max(0, Math.Min(n - mm, jl - ((mm - 2) >> 1)));
            }
            public abstract double rawinterp(int jlo, double x);
        }

        // NumericalRecipes3.Poly_Interp
        public static Tuple<double, double> InterporatePolynomial(IList<double> xx, IList<double> yy, int order, double x)
        {
            double dy = 0.0;
            var c = yy.Take(order).ToArray();
            var d = c.CloneX();
            var ns = xx.MinIndex(v => Math.Abs(v - x));
            var y = yy[ns--];
            for (int m = 1; m < order; m++)
            {
                for (int i = 0; i < order - m; i++)
                {
                    var ho = xx[i] - x;
                    var hp = xx[i + m] - x;
                    var w = c[i + 1] - d[i];
                    if (ho - hp == 0.0) { ThrowException.InvalidOperation("InterporationPoly: divide by 0"); return Tuple.Create(double.NaN, double.NaN); }
                    var den = w / (ho - hp);
                    d[i] = hp * den;
                    c[i] = ho * den;
                }
                dy = 2 * (ns + 1) < (order - m) ? c[ns + 1] : d[ns--];
                y += dy;
            }
            return Tuple.Create(y, dy);
        }

        #endregion

        #region Quadrature

        // NumericalRecipes3.Quadrature
        // NumericalRecipes3.Trapzd
        static Func<double> QuadratureTrapezoid(Func<double, double> function, double start, double end)
        {
            int n = 0;
            double area = 0.0;
            Func<double> next = () =>
            {
                n++;
                if (n == 1)
                {
                    area = 0.5 * (end - start) * (function(start) + function(end));
                }
                else
                {
                    var m = 1 << (n - 2);
                    var w = (end - start) / m;
                    var sum = 0.0;
                    for (int i = 0; i < m; i++)
                        sum += function(start + w * (i + 0.5));
                    area = (area + w * sum) * 0.5;
                }
                return area;
            };
            return next;
        }

        // NumericalRecipes3.qtrap
        public static double IntegrateTrapezoid(Func<double, double> function, double start, double end, double tolerance = 1e-10)
        {
            const int JMAX = 20;
            var quadrature = QuadratureTrapezoid(function, start, end);
            double area = 0.0;
            for (int i = 0; ; i++)
            {
                if (i == JMAX) { ThrowException.WriteLine("IntegrateTrapezoid: too many steps"); break; }
                var old = area;
                area = quadrature();
                if (i > 5)
                    if (Math.Abs(area - old) < tolerance * Math.Abs(old) || (area == 0.0 && old == 0.0)) break;
            }
            return area;
        }

        // NumericalRecipes3.qsimp
        public static double IntegrateSimpson(Func<double, double> function, double start, double end, double tolerance = 1e-10)
        {
            const int JMAX = 20;
            var quadrature = QuadratureTrapezoid(function, start, end);
            double area = 0.0, orig = 0.0;
            for (int i = 0; ; i++)
            {
                if (i == JMAX) { ThrowException.WriteLine("IntegrateSimpson: too many steps"); break; }
                var oldo = orig;
                var olda = area;
                orig = quadrature();
                area = (4.0 * orig - oldo) / 3.0;
                if (i > 5)
                    if (Math.Abs(area - olda) < tolerance * Math.Abs(olda) || (area == 0.0 && olda == 0.0)) break;
            }
            return area;
        }

        // NumericalRecipes3.Midpnt
        static Func<double> QuadratureMidpoint(Func<double, double> function, double start, double end)
        {
            int n = 0;
            double area = 0.0;
            Func<double> next = () =>
            {
                n++;
                if (n == 1)
                {
                    area = (end - start) * function(0.5 * (start + end));
                }
                else
                {
                    int m = (int)Math.Pow(3, n - 2);
                    var w = (end - start) / m;
                    var sum = 0.0;
                    for (int i = 0; i < m; i++)
                    {
                        sum += function(start + w * (i + 1 / 6.0));
                        sum += function(start + w * (i + 5 / 6.0));
                    }
                    area = (area + w * sum) / 3.0;
                }
                return area;
            };
            return next;
        }
        static Func<double[]> QuadratureMidpoint(Func<double, double[]> function, double start, double end)
        {
            int n = 0;
            double[] area = null;
            Func<double[]> next = () =>
            {
                n++;
                if (n == 1)
                {
                    area = Mt.Mul(function(0.5 * (start + end)), end - start);
                }
                else
                {
                    int m = (int)Math.Pow(3, n - 2);
                    var w = (end - start) / m;
                    var sum = new double[area.Length];
                    for (int i = 0; i < m; i++)
                    {
                        sum.LetAdd(function(start + w * (i + 1 / 6.0)));
                        sum.LetAdd(function(start + w * (i + 5 / 6.0)));
                    }
                    area.LetAddMul(sum, w).LetDiv(3);
                }
                return area;
            };
            return next;
        }

        // NumericalRecipes3.Midinf
        static Func<double> QuadratureMidinf(Func<double, double> function, double start, double end)
        {
            Func<double, double> func = x => function(1 / x) / (x * x);
            return QuadratureMidpoint(func, 1 / end, 1 / start);
        }
        // NumericalRecipes3.Midsql
        static Func<double> QuadratureMidsql(Func<double, double> function, double start, double end)
        {
            Func<double, double> func = x => 2.0 * x * function(start + x * x);
            return QuadratureMidpoint(func, 0, Math.Sqrt(end - start));
        }
        // NumericalRecipes3.Midsqu
        static Func<double> QuadratureMidsqu(Func<double, double> function, double start, double end)
        {
            Func<double, double> func = x => 2.0 * x * function(end - x * x);
            return QuadratureMidpoint(func, 0, Math.Sqrt(end - start));
        }
        // NumericalRecipes3.Midexp
        static Func<double> QuadratureMidexp(Func<double, double> function, double start, double end)
        {
            Func<double, double> func = x => function(-Math.Log(x)) / x;
            return QuadratureMidpoint(func, 0, Math.Exp(-start));
        }

        // NumericalRecipes3.qromb
        // NumericalRecipes3.qromo
        static double IntegrateRomberg(Func<double> quadrature, double tolerance, int div, int JMAX)
        {
            const int Order = 5;
            var xxxx = new List<double>(JMAX);
            var area = new List<double>(JMAX);
            var estimated = 0.0;
            for (int i = 0; ; i++)
            {
                if (i == JMAX) { ThrowException.WriteLine("IntegrateRomberg: too many steps"); break; }
                xxxx.Insert(0, 1 / Math.Pow(div, i));
                area.Insert(0, quadrature());
                if (i + 1 >= Order)
                {
                    var r = InterporatePolynomial(xxxx, area, Order, 0.0);
                    estimated = r.Item1;
                    var error = r.Item2;
                    if (Math.Abs(error) <= Math.Abs(estimated) * tolerance) break;
                }
            }
            return estimated;
        }
        public static double IntegrateRombergTrapezoid(Func<double, double> function, double start, double end, double tolerance = 1e-10)
        {
            var quadrature = QuadratureTrapezoid(function, start, end);
            return IntegrateRomberg(quadrature, tolerance, 4, 20);
        }
        public static double IntegrateRombergMidpoint(Func<double, double> function, double start, double end, double tolerance = 3e-9)
        {
            var quadrature = QuadratureMidpoint(function, start, end);
            return IntegrateRomberg(quadrature, tolerance, 9, 14);
        }

        // vector version
        static double[] IntegrateRomberg(Func<double[]> quadrature, double tolerance, int div, int JMAX)
        {
            const int Order = 5;
            var xxxx = new List<double>(JMAX);
            var area = new List<double[]>(JMAX);
            double[] estimated = null;
            for (int i = 0; ; i++)
            {
                if (i == JMAX) { ThrowException.WriteLine("IntegrateRomberg: too many steps"); break; }
                xxxx.Insert(0, 1 / Math.Pow(div, i));
                area.Insert(0, quadrature().CloneX());
                int D = area[0].Length;
                if (i + 1 >= Order)
                {
                    var r = New.Array(D, d => InterporatePolynomial(xxxx, area.Select(v => v[d]).ToArray(), Order, 0.0));
                    estimated = r.Select(v => v.Item1).ToArray();
                    var error = r.Select(v => v.Item2).ToArray();
                    if (Ex.Range(D).All(d => Math.Abs(error[d]) <= Math.Abs(estimated[d]) * tolerance)) break;
                }
            }
            return estimated;
        }
        public static double[] IntegrateRombergMidpoint(Func<double, double[]> function, double start, double end, double tolerance = 3e-9)
        {
            var quadrature = QuadratureMidpoint(function, start, end);
            return IntegrateRomberg(quadrature, tolerance, 9, 14);
        }
        #endregion

        #region Signal processing functions

        #region DataWindow
        public enum DataWindowType { Box, Hanning, Hamming, Blackman, Parzen, Welch, NormalDistribution = -1 }
        public static double[] GetDataWindow(DataWindowType type, int size)
        {
            double[] Table = new double[size];
            Table = new double[size];

            double h = Mt.PI2 / (size - 1);
            switch (type)
            {
                case DataWindowType.Box:
                    for (int i = size; --i >= 0; ) Table[i] = 1;
                    break;
                case DataWindowType.Hanning:
                    for (int i = size; --i >= 0; ) Table[i] = 0.50 - 0.50 * Math.Cos(h * i);
                    break;
                case DataWindowType.Hamming:
                    for (int i = size; --i >= 0; ) Table[i] = 0.54 - 0.46 * Math.Cos(h * i);
                    break;
                case DataWindowType.Blackman:
                    for (int i = size; --i >= 0; ) Table[i] = 0.42 - 0.50 * Math.Cos(h * i) + 0.08 * Math.Cos(2 * h * i);
                    break;
                case DataWindowType.Parzen:
                    for (int i = size; --i >= 0; ) Table[i] = 1.0 - Math.Abs((i * 2 - (size - 1)) / (double)(size + 1));
                    break;
                case DataWindowType.Welch:
                    for (int i = size; --i >= 0; ) Table[i] = 1.0 - Mt.Sq((i * 2 - (size - 1)) / (double)(size + 1));
                    break;
            }
            //double c = Math.Sqrt(n / Table.Sum(x => Sq(x)));
            //for (int i = n; --i >= 0; ) Table[i] *= c;
            return Table;
        }

        public static void Windowing(double[] data, DataWindowType type)
        {
            double[] window = GetDataWindow(type, data.Length);
            for (int i = data.Length; --i >= 0; ) data[i] *= window[i];
        }
        #endregion

        #region Discrete Fourier Transform
        static void FftFactor2_IdealCode0(complex[] data, int isign, int step)
        {
            int n = data.Length;
            var d = Mt.Phase(isign * Math.PI / step);
            var w = complex.One;
            for (int k = 0; k < step; k++)
            {
                for (int i = k; i < n; i += step * 2)  //stepping access, low calculation
                {
                    int j = i + step;
                    complex t = data[j] * w;  //w = Mt.Phase(isign * Math.PI / step * k);
                    data[j] = data[i] - t;
                    data[i] += t;
                }
                w *= d;
            }
        }
        unsafe static void FftFactor2_TestCode0(complex* datacomplex, int n, int isign, int step)
        {
            double* data = (double*)datacomplex;
            double theta = isign * Math.PI / step;
            double ddxr = Math.Cos(theta);
            double ddxi = Math.Sin(theta);
            double dxr = 1.0;
            double dxi = 0.0;
            for (int k = 0; k < step; k++)
            {
                for (int i = k; i < n; i += step * 2)  //stepping access, low calculation
                {
                    int j = i + step;
                    double tr = data[j * 2 + 0] * dxr - data[j * 2 + 1] * dxi;
                    double ti = data[j * 2 + 0] * dxi + data[j * 2 + 1] * dxr;
                    data[j * 2 + 0] = data[i * 2 + 0] - tr;
                    data[j * 2 + 1] = data[i * 2 + 1] - ti;
                    data[i * 2 + 0] += tr;
                    data[i * 2 + 1] += ti;
                }
                double _r = dxr;
                dxr = _r * ddxr - dxi * ddxi;
                dxi = _r * ddxi + dxi * ddxr;
            }
        }
        static void FftFactor2_IdealCode1(complex[] data, int isign, int step)
        {
            int n = data.Length;
            var d = Mt.Phase(isign * Math.PI / step);
            for (int h = 0; h < n; h += step * 2)
            {
                var w = complex.One;
                for (int k = 0; k < step; k++)  //sequential access, high calculation
                {
                    int i = h + k;
                    int j = i + step;
                    complex t = data[j] * w;  //w = Mt.Phase((Mt.PI2 / n) * (isign * (n / (step * 2))) * k);
                    data[j] = data[i] - t;
                    data[i] += t;
                    w *= d;
                }
            }
        }
        unsafe static void FftFactor2(complex* datacomplex, int n, int isign, int step)
        {
            double* data = (double*)datacomplex;
            double theta = isign * Math.PI / step;
            double ddxr = Math.Cos(theta);
            double ddxi = Math.Sin(theta);
            for (int h = 0; h < n; h += step * 2)
            {
                double dxr = 1.0;
                double dxi = 0.0;
                for (int k = 0; k < step; k++)  //sequential access, high calculation
                {
                    int i = h + k;
                    int j = i + step;
                    double tr = data[j * 2 + 0] * dxr - data[j * 2 + 1] * dxi;
                    double ti = data[j * 2 + 0] * dxi + data[j * 2 + 1] * dxr;
                    data[j * 2 + 0] = data[i * 2 + 0] - tr;
                    data[j * 2 + 1] = data[i * 2 + 1] - ti;
                    data[i * 2 + 0] += tr;
                    data[i * 2 + 1] += ti;
                    double _r = dxr;
                    dxr = _r * ddxr - dxi * ddxi;
                    dxi = _r * ddxi + dxi * ddxr;
                }
            }
        }

        static void FftFactorN_IdealCode1(complex[] data, int isign, int step, int factor)
        {
            int n = data.Length;
            var temp = new complex[factor];
            var ddx = Mt.Phase(isign * Mt.PI2 / (step * factor));
            var dy = Mt.Phase(isign * Mt.PI2 / factor);
            for (int h = 0; h < n; h += step * factor)
            {
                var dx = complex.One;
                for (int s = 0; s < step; s++)
                {
                    for (int ff = 0; ff < factor; ff++) temp[ff] = data[h + s];
                    var x = dx;
                    var y = dy;
                    for (int f = 1; f < factor; f++)
                    {
                        complex z = data[h + s + f * step] * x;  //x = s * f
                        for (int ff = 0; ff < factor; ff++)
                        {
                            temp[ff] += z;  //z = data[h + s + f * step] * Mt.Phase(isign * Mt.PI2 / (step * factor) * (s + step * ff) * f);
                            z *= y;  //y = step * f
                        }
                        y *= dy;
                        x *= dx;
                    }
                    for (int f = 0; f < factor; f++) data[h + s + f * step] = temp[f];
                    dx *= ddx;
                }
            }
        }
        unsafe static void FftFactorN_TestCode1(complex* datacomplex, int n, int isign, int step, int factor)
        {
            fixed (complex* tempcomplex = new complex[factor])
            {
                double* data = (double*)datacomplex;
                double* temp = (double*)tempcomplex;
                double theta = isign * Mt.PI2 / (step * factor);
                double ddxr = Math.Cos(theta);
                double ddxi = Math.Sin(theta);
                double rho = isign * Mt.PI2 / factor;
                double dyr = Math.Cos(rho);
                double dyi = Math.Sin(rho);
                for (int h = 0; h < n; h += step * factor)
                {
                    double dxr = 1.0;
                    double dxi = 0.0;
                    for (int s = 0; s < step; s++)
                    {
                        for (int ff = 0; ff < factor; ff++) Us.Let(&tempcomplex[ff], &datacomplex[h + s]);
                        double xr = dxr;
                        double xi = dxi;
                        double yr = dyr;
                        double yi = dyi;
                        double _r;
                        for (int f = 1; f < factor; f++)
                        {
                            int i = h + s + f * step;
                            double zr = data[i * 2 + 0] * xr - data[i * 2 + 1] * xi;
                            double zi = data[i * 2 + 0] * xi + data[i * 2 + 1] * xr;
                            for (int ff = 0; ff < factor; ff++)
                            {
                                temp[ff * 2 + 0] += zr;
                                temp[ff * 2 + 1] += zi;
                                _r = zr;
                                zr = _r * yr - zi * yi;
                                zi = _r * yi + zi * yr;
                            }
                            _r = yr;
                            yr = _r * dyr - yi * dyi;
                            yi = _r * dyi + yi * dyr;
                            _r = xr;
                            xr = _r * dxr - xi * dxi;
                            xi = _r * dxi + xi * dxr;
                        }
                        for (int f = 0; f < factor; f++) Us.Let(&datacomplex[h + s + f * step], &tempcomplex[f]);
                        _r = dxr;
                        dxr = _r * ddxr - dxi * ddxi;
                        dxi = _r * ddxi + dxi * ddxr;
                    }
                }
            }
        }
        static void FftFactorN_IdealCode2(complex[] data, int isign, int step, int factor)
        {
            int n = data.Length;
            if (step > 1)
            {
                var ddx = Mt.Phase(isign * Mt.PI2 / (step * factor));
                for (int h = 0; h < n; h += step * factor)
                {
                    var dx = complex.One;
                    for (int f = 1; f < factor; f++)
                    {
                        dx *= ddx;
                        var x = complex.One;
                        for (int s = 1; s < step; s++)
                        {
                            x *= dx;
                            data[h + s + f * step] *= x;  //x = Mt.Phase(isign * Mt.PI2 / (step * factor) * (s * f))
                        }
                    }
                }
            }
            var temp = new complex[factor];
            var dy = Mt.Phase(isign * Mt.PI2 / factor);
            for (int h = 0; h < n; h += step * factor)
            {
                for (int s = 0; s < step; s++)
                {
                    for (int f = 0; f < factor; f++) temp[f] = data[h + s];
                    var y = complex.One;
                    for (int f = 1; f < factor; f++)
                    {
                        y *= dy;
                        var z = data[h + s + f * step];
                        for (int ff = 0; ff < factor; ff++)
                        {
                            temp[ff] += z;
                            z *= y;
                        }
                    }
                    for (int f = 0; f < factor; f++) data[h + s + f * step] = temp[f];
                }
            }
        }
        unsafe static void FftFactorN(complex* datacomplex, int n, int isign, int step, int factor)
        {
            double* data = (double*)datacomplex;
            double _r, _i;
            if (step > 1)
            {
                double theta = isign * Mt.PI2 / (step * factor);
                double ddxr = Math.Cos(theta);
                double ddxi = Math.Sin(theta);
                for (int h = 0; h < n; h += step * factor)
                {
                    double dxr = 1.0;
                    double dxi = 0.0;
                    for (int f = 1; f < factor; f++)
                    {
                        _r = dxr;
                        dxr = _r * ddxr - dxi * ddxi;
                        dxi = _r * ddxi + dxi * ddxr;
                        double xr = 1.0;
                        double xi = 0.0;
                        for (int k = 1; k < step; k++)
                        {
                            _r = xr;
                            xr = _r * dxr - xi * dxi;
                            xi = _r * dxi + xi * dxr;
                            int i = h + k + f * step;
                            _r = data[i * 2 + 0];
                            _i = data[i * 2 + 1];
                            data[i * 2 + 0] = _r * xr - _i * xi;
                            data[i * 2 + 1] = _r * xi + _i * xr;
                        }
                    }
                }
            }
            fixed (complex* tempcomplex = new complex[factor])
            {
                double* temp = (double*)tempcomplex;
                double rho = isign * Mt.PI2 / factor;
                double ddyr = Math.Cos(rho);
                double ddyi = Math.Sin(rho);
                for (int h = 0; h < n; h += step * factor)
                {
                    for (int s = 0; s < step; s++)
                    {
                        for (int f = 0; f < factor; f++) Us.Let(&tempcomplex[f], &datacomplex[h + s]);
                        double dyr = 1.0;
                        double dyi = 0.0;
                        for (int f = 1; f < factor; f++)
                        {
                            _r = dyr;
                            dyr = _r * ddyr - dyi * ddyi;
                            dyi = _r * ddyi + dyi * ddyr;
                            int i = h + s + f * step;
                            double zr = data[i * 2 + 0];
                            double zi = data[i * 2 + 1];
                            temp[0] += zr;
                            temp[1] += zi;
                            for (int ff = 1; ff < factor; ff++)
                            {
                                _r = zr;
                                temp[ff * 2 + 0] += zr = _r * dyr - zi * dyi;
                                temp[ff * 2 + 1] += zi = _r * dyi + zi * dyr;
                            }
                        }
                        for (int f = 0; f < factor; f++) Us.Let(&datacomplex[h + s + f * step], &tempcomplex[f]);
                    }
                }
            }
        }

        public static Func<int, int[]> FactorIntegerExpandedReverse = New.Cache(Mt.FactorIntegerExpanded);

        static void FFT_IdealCode(complex[] data, bool forward)
        {
            if (data == null) ThrowException.Argument("FFT: data");
            int n = data.Length;
            if (n < 0) ThrowException.Argument("FFT: n");
            var isign = forward ? -1 : +1;  //NumericalRecipesはisignの意味が逆なので注意
            if ((n & (n - 1)) == 0)
            {
                for (int j = 0, i = 0; i < n; i++)
                {
                    if (i < j) Mt.Swap(ref data[i], ref data[j]);
                    for (int m = n; j < m; j ^= (m >>= 1)) ;
                }
                for (int step = 1; step < n; step *= 2)
                    FftFactor2_IdealCode1(data, isign, step);
            }
            else
            {
                var factors = FactorIntegerExpandedReverse(n);
                {
                    var dim = factors.Length;
                    var mm = new int[dim];
                    var ii = new int[dim];
                    for (int nn = 1, d = dim; --d >= 0; ) { mm[d] = nn; nn *= factors[d]; }
                    int j = 0;
                    var temp = new complex[n];
                    for (int i = 0; i < n; i++)
                    {
                        temp[i] = data[j];
                        for (int d = 0; d < dim; d++)
                        {
                            if (++ii[d] < factors[d]) { j += mm[d]; break; }
                            ii[d] = 0; j -= mm[d] * (factors[d] - 1);
                        }
                    }
                    Array.Copy(temp, data, n);
                }
                for (int d = 0, step = 1; step < n; d++)
                {
                    var factor = factors[d];
                    switch (factor)
                    {
                        case 2: FftFactor2_IdealCode1(data, isign, step); break;
                        default: FftFactorN_IdealCode1(data, isign, step, factor); break;
                    }
                    step *= factor;
                }
            }
        }
        public static complex[] InverseDiscreteFourierTransform_IdealCode(complex[] data) { return DiscreteFourierTransform_IdealCode(data, false); }
        public static complex[] DiscreteFourierTransform_IdealCode(complex[] data, bool forward = true)
        {
            var result = data.CloneX();
            FFT_IdealCode(result, forward);
            result.LetDiv(Math.Sqrt(result.Length));  //unitary transform
            return result;
        }

        //必ず通過
        unsafe static void FFT(complex* data, int n, bool forward)
        {
            if (data == null) ThrowException.Argument("FFT: data");
            if (n < 0) ThrowException.Argument("FFT: n");
            var isign = forward ? -1 : +1;  //NumericalRecipesはisignの意味が逆なので注意
            if ((n & (n - 1)) == 0)
            {
                for (int j = 0, i = 0; i < n; i++)
                {
                    if (i < j) Us.Swap(&data[i], &data[j]);
                    for (int m = n; j < m; j ^= (m >>= 1)) ;
                }
                for (int step = 1; step < n; step *= 2)
                    FftFactor2(data, n, isign, step);
            }
            else
            {
                var factors = FactorIntegerExpandedReverse(n);
                var tempArray = new complex[n];
                fixed (complex* temp = tempArray)
                {
                    var dim = factors.Length;
                    var mm = new int[dim];
                    var ii = new int[dim];
                    for (int nn = 1, d = dim; --d >= 0; ) { mm[d] = nn; nn *= factors[d]; }
                    for (int j = 0, i = 0; i < n; i++)
                    {
                        temp[i] = data[j];
                        for (int d = 0; d < dim; d++)
                        {
                            if (++ii[d] < factors[d]) { j += mm[d]; break; }
                            ii[d] = 0; j -= mm[d] * (factors[d] - 1);
                        }
                    }
                    Us.Let(data, temp, n);
                }
                for (int d = 0, step = 1; step < n; d++)
                {
                    var factor = factors[d];
                    switch (factor)
                    {
                        case 2: FftFactor2(data, n, isign, step); break;
                        default: FftFactorN(data, n, isign, step, factor); break;
                    }
                    step *= factor;
                }
            }
        }
        unsafe static void FFT(complex* data, int[] nn, bool forward)
        {
            int ntot = nn.Product();
            if (nn.Length == 1)
            {
                FFT(data, ntot, forward);
            }
            else
            {
                for (int d = nn.Length, step = 1; --d >= 0; )
                {
                    int n = nn[d];
                    var tempArray = new complex[n];
                    fixed (complex* temp = tempArray)
                    {
                        for (int i2 = 0; i2 < ntot; i2 += n * step)
                        {
                            for (int i1 = 0; i1 < step; i1++)
                            {
                                Us.LetStartStep(temp, data, i1 + i2, step, n);
                                FFT(temp, n, forward);
                                Us.LetStartStep(data, i1 + i2, step, temp, n);
                            }
                        }
                    }
                    step *= n;
                }
            }
            Us.Mul(data, data, 1 / Math.Sqrt(ntot), ntot);  //unitary transform
        }

        public static complex[] InverseDiscreteFourierTransform(complex[] data) { return DiscreteFourierTransform(data, false); }
        public static complex[] DiscreteFourierTransform(complex[] data, bool forward = true)
        {
            var result = data.CloneX();
            unsafe { fixed (complex* r = result) FFT(r, data.GetLengths(), forward); }
            return result;
        }
        public static complex[,] InverseDiscreteFourierTransform(complex[,] data) { return DiscreteFourierTransform(data, false); }
        public static complex[,] DiscreteFourierTransform(complex[,] data, bool forward = true)
        {
            var result = data.CloneX();
            unsafe { fixed (complex* r = result) FFT(r, data.GetLengths(), forward); }
            return result;
        }
        public static complex[, ,] InverseDiscreteFourierTransform(complex[, ,] data) { return DiscreteFourierTransform(data, false); }
        public static complex[, ,] DiscreteFourierTransform(complex[, ,] data, bool forward = true)
        {
            var result = data.CloneX();
            unsafe { fixed (complex* r = result) FFT(r, data.GetLengths(), forward); }
            return result;
        }

        public static complex[] InverseDiscreteFourierTransform(double[] data) { return DiscreteFourierTransform(data, false); }
        public static complex[] DiscreteFourierTransform(double[] data, bool forward = true)
        {
            return DiscreteFourierTransform(data.Tocomplex(), forward);
        }

        public static complex[] InverseDiscreteFourierTransform_DefinitionCode(complex[] data) { return DiscreteFourierTransform_DefinitionCode(data, false); }
        public static complex[] DiscreteFourierTransform_DefinitionCode(complex[] data, bool forward = true)
        {
            int n = data.Length;
            int isign = forward ? -1 : +1;
            var s = 1 / Math.Sqrt(data.Length);
            return New.Array(n, i => Mt.Sum(n, j => data[j] * Mt.Phase(isign * i * j, n)) * s);
        }

        //NumericalRecipes3.realft
        //bug-fixed
        unsafe static void RealFFT(double* data, int n, bool forward)
        {
            if ((n & 3) != 0) ThrowException.Argument("RealFFT: n should be multiple of 4");
            int n2 = n >> 1;
            if (forward)
            {
                FFT((complex*)data, n2, forward);
            }
            double theta = (forward ? +1 : -1) * Mt.PI2 / n;
            double dr = Math.Cos(theta);
            double di = Math.Sin(theta);
            double wr = 1.0;
            double wi = 0.0;
            for (int i = (n >> 2); i > 0; i--)
            {
                int j = n2 - i;
                double h1r = data[i * 2 + 0] + data[j * 2 + 0];
                double h1i = data[i * 2 + 1] - data[j * 2 + 1];
                double h2r = data[i * 2 + 0] - data[j * 2 + 0];
                double h2i = data[i * 2 + 1] + data[j * 2 + 1];
                double tr = h2r * wr - h2i * wi;
                double ti = h2r * wi + h2i * wr;
                data[i * 2 + 0] = +0.5 * (h1r - tr);
                data[i * 2 + 1] = +0.5 * (h1i - ti);
                data[j * 2 + 0] = +0.5 * (h1r + tr);
                data[j * 2 + 1] = -0.5 * (h1i + ti);
                double _r = wr;
                wr = _r * dr - wi * di;
                wi = _r * di + wi * dr;
            }
            double d0 = data[0];
            data[0] = d0 + data[1];
            data[1] = d0 - data[1];
            if (!forward)
            {
                data[0] *= 0.5;
                data[1] *= 0.5;
                FFT((complex*)data, n2, forward);
            }
        }
        public unsafe static complex[] RealDiscreteFourierTransform(double[] data)
        {
            int n = data.Length;
            var result = new complex[n / 2 + 1];
            fixed (complex* resultcomplex = result)
            fixed (double* d = data)
            {
                double* r = (double*)resultcomplex;
                Us.Mul(r, d, 1 / Math.Sqrt(n), n);  //unitary transform
                RealFFT(r, n, true);
                r[n] = r[1];
                r[1] = 0;
            }
            return result;
        }
        public unsafe static double[] RealInverseDiscreteFourierTransform(complex[] data)
        {
            int n = (data.Length - 1) * 2;
            var result = new double[n];
            fixed (double* r = result)
            fixed (complex* datacomplex = data)
            {
                double* d = (double*)datacomplex;
                double c = 2 / Math.Sqrt(n);  //unitary transform
                Us.Mul(r, d, c, n);
                r[1] = d[n] * c;
                RealFFT(r, n, false);
            }
            return result;
        }

        public static double[] PowerSpectrumFFT(complex[] data)
        {
            var fft = DiscreteFourierTransform(data);
            var result = New.Array(fft.Length / 2 + 1, i => fft[i].SqAbs() * 2);
            result[0] /= 2;
            if ((fft.Length & 1) == 0) result[result.Length - 1] /= 2;
            return result;
        }

        #endregion

        #region Convolution
        public enum ConvolutionType { Full, Same, Valid };
        public static int ConvolutionSize(int length0, int length1, ConvolutionType type)
        {
            switch (type)
            {
                case ConvolutionType.Full: return Math.Max(0, length0 + length1 - 1);
                case ConvolutionType.Same: return length0;
                case ConvolutionType.Valid: return Math.Max(0, length0 - Math.Max(0, length1 - 1));
            }
            return ThrowException.Argument<int>("ConvolutionSize: type is unknown");
        }

        public static complex[] Convolution_DefinitionCode(complex[] data0, complex[] data1, ConvolutionType type = ConvolutionType.Full)
        {
            var n0 = data0.Length;
            var n1 = data1.Length;
            var nR = ConvolutionSize(n0, n1, type);
            var dataR = new complex[nR];
            for (int i0 = 0; i0 < n0; i0++)
                for (int i1 = 0; i1 < n1; i1++)
                {
                    var j = i0 + i1 - (n0 + n1 - nR) / 2;
                    if (j >= 0 && j < nR) dataR[j] += data0[i0] * data1[i1];
                }
            return dataR;
        }

        static IEnumerable<Int2> CopyToEnlargeIterator(int[] nn0, int[] nn1)
        {
            var i = default(Int2);
            var dim = nn0.Length - 1;
            if (dim == 0) { yield return i; yield break; }
            var indx = new int[dim];
            var skip = new int[dim];
            var nnn = 1;
            for (int d = dim; --d >= 0; nnn *= nn1[d]) skip[d] = nnn * (nn1[d] - nn0[d]);
            for (; i.v1 < nnn; i.v0++, i.v1++)
            {
                yield return i;
                for (int d = dim; --d >= 0 && ++indx[d] == nn0[d]; ) { indx[d] = 0; i.v1 += skip[d]; }
            }
        }
        unsafe static void CopyToEnlarge(complex* data0, int[] nn0, complex* data1, int[] nn1)
        {
            int n0 = nn0.Last();
            int n1 = nn1.Last();
            foreach (var i in CopyToEnlargeIterator(nn0, nn1))
                Us.Let(data1 + i.v1 * n1, data0 + i.v0 * n0, n0);
        }
        static IEnumerable<Int2> CopyToEnlargeShrinkIterator(int[] nn0, int[] nn1, int[] oo1)
        {
            var i = default(Int2);
            var dim = nn0.Length - 1;
            if (dim == 0) { yield return i; yield break; }
            var indx = new int[dim];
            var skip = new int[dim];
            var nnn = 1;
            for (int d = dim; --d >= 0; nnn *= nn0[d]) { i.v0 += nnn * oo1[d]; skip[d] = nnn * (nn0[d] - nn1[d] + oo1[d]); }
            for (; i.v0 < nnn; i.v0++, i.v1++)
            {
                yield return i;
                for (int d = dim; --d >= 0 && ++indx[d] == nn1[d]; ) { indx[d] = 0; i.v0 += skip[d]; }
            }
        }
        unsafe static void CopyToShrink(complex* data0, int[] nn0, complex* data1, int[] nn1, int[] oo1)
        {
            int n0 = nn0.Last();
            int n1 = nn1.Last();
            int o0 = oo1.Last();
            foreach (var i in CopyToEnlargeShrinkIterator(nn0, nn1, oo1))
                Us.Let(data1 + i.v1 * n1, data0 + (i.v0 * n0 + o0), n1);
        }

        unsafe static Array Convolution(Array data0, Array data1, ConvolutionType type, bool filter)
        {
            if (data0.GetType().GetElementType() != typeof(complex)) return null;
            var nn0 = data0.GetLengths();
            var nn1 = data1.GetLengths();
            var nnR = nn0.ZipTo(nn1, (i0, i1) => ConvolutionSize(i0, i1, type));
            complex* d0 = null, d1 = null, dR = null;
            Action action = () =>
            {
                if (data0.Length == 0 || data1.Length == 0) return;
                var nnT = nn0.ZipTo(nn1, (i0, i1) => 1 << Mt.Log2Ceiling(Math.Max(1, i0 + i1 - 1)));
                var ooR = nn0.ZipTo(nn1, nnR, (i0, i1, i2) => (i0 + i1 - i2) / 2);
                var nT = nnT.Product();

                fixed (complex* t0 = new complex[nT], t1 = new complex[nT], tF = filter ? new complex[data1.Length] : null)
                {
                    if (filter) { Us.Rev(tF, d1, data1.Length); d1 = tF; }
                    CopyToEnlarge(d0, nn0, t0, nnT); Nm.FFT(t0, nnT, true);
                    CopyToEnlarge(d1, nn1, t1, nnT); Nm.FFT(t1, nnT, true);
                    Us.Mul(t0, t0, t1, nT);
                    Us.Mul(t0, t0, Math.Sqrt(nT), nT);
                    Nm.FFT(t0, nnT, false); CopyToShrink(t0, nnT, dR, nnR, ooR);
                }
            };
            var dataR = Array.CreateInstance(typeof(complex), nnR);
            switch (data0.Rank)
            {
                case 1: fixed (complex* _d0 = (complex[])data0, _d1 = (complex[])data1, _dR = (complex[])dataR) { d0 = _d0; d1 = _d1; dR = _dR; action(); } break;
                case 2: fixed (complex* _d0 = (complex[,])data0, _d1 = (complex[,])data1, _dR = (complex[,])dataR) { d0 = _d0; d1 = _d1; dR = _dR; action(); } break;
                case 3: fixed (complex* _d0 = (complex[, ,])data0, _d1 = (complex[, ,])data1, _dR = (complex[, ,])dataR) { d0 = _d0; d1 = _d1; dR = _dR; action(); } break;
            }
            return dataR;
        }
        public unsafe static Array Convolution(Array data0, Array data1, ConvolutionType type = ConvolutionType.Full) { return Convolution(data0, data1, type, false); }
        public static complex[] Convolution(complex[] data0, complex[] data1, ConvolutionType type = ConvolutionType.Full) { return (complex[])Convolution((Array)data0, (Array)data1, type); }
        public static complex[,] Convolution(complex[,] data0, complex[,] data1, ConvolutionType type = ConvolutionType.Full) { return (complex[,])Convolution((Array)data0, (Array)data1, type); }
        public static complex[, ,] Convolution(complex[, ,] data0, complex[, ,] data1, ConvolutionType type = ConvolutionType.Full) { return (complex[, ,])Convolution((Array)data0, (Array)data1, type); }
        public static double[] Convolution(double[] data0, double[] data1, ConvolutionType type = ConvolutionType.Full) { return Convolution(data0.Tocomplex(), data1.Tocomplex(), type).Real(); }
        public static double[,] Convolution(double[,] data0, double[,] data1, ConvolutionType type = ConvolutionType.Full) { return Convolution(data0.Tocomplex(), data1.Tocomplex(), type).Real(); }
        public static double[, ,] Convolution(double[, ,] data0, double[, ,] data1, ConvolutionType type = ConvolutionType.Full) { return Convolution(data0.Tocomplex(), data1.Tocomplex(), type).Real(); }

        public unsafe static Array Filter(Array data, Array filter, ConvolutionType type = ConvolutionType.Same) { return Convolution(data, filter, type, true); }
        public static complex[] Filter(complex[] data, complex[] filter, ConvolutionType type = ConvolutionType.Same) { return (complex[])Filter((Array)data, (Array)filter, type); }
        public static complex[,] Filter(complex[,] data, complex[,] filter, ConvolutionType type = ConvolutionType.Same) { return (complex[,])Filter((Array)data, (Array)filter, type); }
        public static complex[, ,] Filter(complex[, ,] data, complex[, ,] filter, ConvolutionType type = ConvolutionType.Same) { return (complex[, ,])Filter((Array)data, (Array)filter, type); }
        public static double[] Filter(double[] data, double[] filter, ConvolutionType type = ConvolutionType.Same) { return Filter(data.Tocomplex(), filter.Tocomplex(), type).Real(); }
        public static double[,] Filter(double[,] data, double[,] filter, ConvolutionType type = ConvolutionType.Same) { return Filter(data.Tocomplex(), filter.Tocomplex(), type).Real(); }
        public static double[, ,] Filter(double[, ,] data, double[, ,] filter, ConvolutionType type = ConvolutionType.Same) { return Filter(data.Tocomplex(), filter.Tocomplex(), type).Real(); }
        #endregion

        #region Discrete Sine and Cosine Transform
        // NumericalRecipes3.dct2
        // modified to meet the definition of FFT/real FFT
        unsafe static void DCT2(double* y, int n, bool forward)
        {
            double theta = Math.PI / n;
            double wpr = Math.Cos(theta);
            double wpi = Math.Sin(theta);
            double wqr = Math.Cos(theta / 2);
            double wqi = Math.Sin(theta / 2);
            double wr = 1.0;
            double wi = 0.0;
            double _r;
            if (forward)
            {
                for (int i = 0; i < n / 2; i++)
                {
                    int j = n - 1 - i;
                    double y1 = (y[i] + y[j]) * 0.5;
                    double y2 = (y[i] - y[j]) * wqi;
                    y[i] = y1 + y2;
                    y[j] = y1 - y2;
                    _r = wqr;
                    wqr = _r * wpr - wqi * wpi;
                    wqi = _r * wpi + wqi * wpr;
                }
                RealFFT(y, n, true);  //複素共役
                for (int i = 2; i < n; i += 2)
                {
                    _r = wr;
                    wr = _r * wpr - wi * wpi;
                    wi = _r * wpi + wi * wpr;
                    double y1 = y[i + 0] * wr + y[i + 1] * wi;  //複素共役
                    double y2 = y[i + 0] * wi - y[i + 1] * wr;  //複素共役
                    y[i + 0] = y1;
                    y[i + 1] = y2;
                }
                double sum = 0.5 * y[1];
                for (int i = n - 1; i > 0; i -= 2)
                {
                    double sum1 = sum;
                    sum += y[i];
                    y[i] = sum1;
                }
            }
            else
            {
                double ytemp = y[n - 1];
                for (int i = n - 1; i > 2; i -= 2)
                    y[i] = y[i - 2] - y[i];
                y[1] = ytemp * 2;
                for (int i = 2; i < n; i += 2)
                {
                    _r = wr;
                    wr = _r * wpr - wi * wpi;
                    wi = _r * wpi + wi * wpr;
                    double y1 = y[i + 0] * wr + y[i + 1] * wi;
                    double y2 = y[i + 1] * wr - y[i + 0] * wi;
                    y[i + 0] = +y1;
                    y[i + 1] = -y2;  //複素共役
                }
                RealFFT(y, n, false);  //複素共役
                for (int i = 0; i < n / 2; i++)
                {
                    int j = n - 1 - i;
                    double y1 = (y[i] + y[j]);
                    double y2 = (y[i] - y[j]) * (0.5 / wqi);
                    y[i] = 0.5 * (y1 + y2);
                    y[j] = 0.5 * (y1 - y2);
                    _r = wqr;
                    wqr = _r * wpr - wqi * wpi;
                    wqi = _r * wpi + wqi * wpr;
                }
            }
        }

        public unsafe static double[] DiscreteCosineTransform(double[] data, bool forward = true)
        {
            int n = data.Length;
            var result = data.CloneX();
            fixed (double* r = result)
            {
                if (!forward) r[0] *= Mt.Sqrt2;
                DCT2(r, n, forward);
                if (forward) r[0] /= Mt.Sqrt2;
            }
            return result.LetMul(Math.Sqrt(2.0 / n));
        }

        public static double[] InverseDiscreteCosineTransform_DefinitionCode(double[] data, bool forward = true) { return DiscreteCosineTransform_DefinitionCode(data, false); }
        public static double[] DiscreteCosineTransform_DefinitionCode(double[] data, bool forward = true)
        {
            int n = data.Length;
            return (forward ?
                New.Array(n, k => Mt.Sum(n, i => Math.Cos(Math.PI / n * k * (i + 0.5)) * data[i]) * (k == 0 ? 1 / Mt.Sqrt2 : 1)) :
                New.Array(n, k => Mt.Sum(n, i => Math.Cos(Math.PI / n * (k + 0.5) * i) * data[i] * (i == 0 ? 1 / Mt.Sqrt2 : 1))))
                .LetMul(Math.Sqrt(2.0 / n));
        }

        public static double[] InverseDiscreteCosineTransformType1_DefinitionCode(double[] data) { return DiscreteCosineTransformType1_DefinitionCode(data, false); }
        public static double[] DiscreteCosineTransformType1_DefinitionCode(double[] data, bool forward = true)
        {
            int n = data.Length;
            double c = Math.Sqrt(2.0 / (n - 1));
            return New.Array(n, k => Mt.Sum(n, i => Math.Cos(Math.PI / (n - 1) * k * i) * data[i] * c * ((i == 0 || i == n - 1) ? 0.5 : 1.0)));
        }
        public static double[] InverseDiscreteCosineTransformType2_DefinitionCode(double[] data) { return DiscreteCosineTransformType2_DefinitionCode(data, false); }
        public static double[] DiscreteCosineTransformType2_DefinitionCode(double[] data, bool forward = true)
        {
            int n = data.Length;
            if (!forward) return DiscreteCosineTransformType3_DefinitionCode(data);
            double c = Math.Sqrt(2.0 / n);
            return New.Array(n, k => Mt.Sum(n, i => Math.Cos(Math.PI / n * k * (i + 0.5)) * data[i] * c));
        }
        public static double[] InverseDiscreteCosineTransformType3_DefinitionCode(double[] data) { return DiscreteCosineTransformType3_DefinitionCode(data, false); }
        public static double[] DiscreteCosineTransformType3_DefinitionCode(double[] data, bool forward = true)
        {
            int n = data.Length;
            if (!forward) return DiscreteCosineTransformType2_DefinitionCode(data);
            double c = Math.Sqrt(2.0 / n);
            return New.Array(n, k => Mt.Sum(n, i => Math.Cos(Math.PI / n * (k + 0.5) * i) * data[i] * (i == 0 ? 0.5 : 1) * c));
        }
        public static double[] InverseDiscreteCosineTransformType4_DefinitionCode(double[] data) { return DiscreteCosineTransformType4_DefinitionCode(data, false); }
        public static double[] DiscreteCosineTransformType4_DefinitionCode(double[] data, bool forward = true)
        {
            int n = data.Length;
            double c = Math.Sqrt(2.0 / n);
            return New.Array(n, k => Mt.Sum(n, i => Math.Cos(Math.PI / n * (k + 0.5) * (i + 0.5)) * data[i] * c));
        }

        unsafe static void DCT_DefinitionCode(double* data, int n, bool forward)
        {
            if (!forward) data[0] /= Mt.Sqrt2;
            var result = forward ?
                New.Array(n, k => { var c = Math.PI / n * k; return Mt.Sum(n, i => Math.Cos(c * (i + 0.5)) * data[i]); }) :
                New.Array(n, k => { var c = Math.PI / n * (k + 0.5); return Mt.Sum(n, i => Math.Cos(c * i) * data[i]); });
            if (forward) result[0] /= Mt.Sqrt2;
            fixed (double* r = result)
                Us.Let(data, r, n);
        }
        unsafe static void DCT_DefinitionCode(double* data, int[] nn, bool forward)
        {
            int ntot = nn.Product();
            if (nn.Length == 1)
            {
                DCT_(data, ntot, forward);
            }
            else
            {
                nn.ForEach((n, d) =>
                {
                    int step = Mt.Product(d, i => nn[i]);
                    var tempArray = new double[n];
                    fixed (double* temp = tempArray)
                    {
                        for (int i2 = 0; i2 < ntot; i2 += n * step)
                        {
                            for (int i1 = 0; i1 < step; i1++)
                            {
                                Us.LetStartStep(temp, data, i1 + i2, step, n);
                                DCT_DefinitionCode(temp, n, forward);
                                Us.LetStartStep(data, i1 + i2, step, temp, n);
                            }
                        }
                    }
                });
            }
            Us.Mul(data, data, Math.Sqrt((1 << nn.Length) / (double)ntot), ntot);  //unitary transform
        }
        public static double[,] InverseDiscreteCosineTransform_DefinitionCode(double[,] data) { return DiscreteCosineTransform_DefinitionCode(data, false); }
        public static double[,] DiscreteCosineTransform_DefinitionCode(double[,] data, bool forward = true)
        {
            var result = data.CloneX(); unsafe { fixed (double* r = result) DCT_DefinitionCode(r, data.GetLengths(), forward); }
            return result;
        }

        unsafe static void DCT_(double* data, int n, bool forward)
        {
            if (!forward) data[0] *= Mt.Sqrt2;
            DCT2(data, n, forward);
            if (forward) data[0] /= Mt.Sqrt2;
        }
        unsafe static void DCT(double* data, int[] nn, bool forward)
        {
            int ntot = nn.Product();
            if (nn.Length == 1)
            {
                DCT_(data, ntot, forward);
            }
            else
            {
                for (int d = nn.Length, step = 1; --d >= 0; )
                {
                    int n = nn[d];
                    var tempArray = new double[n];
                    fixed (double* temp = tempArray)
                    {
                        for (int i2 = 0; i2 < ntot; i2 += n * step)
                        {
                            for (int i1 = 0; i1 < step; i1++)
                            {
                                Us.LetStartStep(temp, data, i1 + i2, step, n);
                                DCT_(temp, n, forward);
                                Us.LetStartStep(data, i1 + i2, step, temp, n);
                            }
                        }
                    }
                    step *= n;
                }
            }
            Us.Mul(data, data, Math.Sqrt((1 << nn.Length) / (double)ntot), ntot);  //unitary transform
        }
        public static double[,] InverseDiscreteCosineTransform(double[,] data) { return DiscreteCosineTransform(data, false); }
        public static double[,] DiscreteCosineTransform(double[,] data, bool forward = true)
        {
            var result = data.CloneX(); unsafe { fixed (double* r = result) DCT(r, data.GetLengths(), forward); }
            return result;
        }

        public static complex[] DiscreteCosineTransform(complex[] data, bool forward = true)
        {
            var real = data.Real(); unsafe { fixed (double* r = real) DCT(r, real.GetLengths(), forward); }
            var imag = data.Imag(); unsafe { fixed (double* r = imag) DCT(r, imag.GetLengths(), forward); }
            return Mt.Tocomplex(real, imag);
        }
        public static complex[,] DiscreteCosineTransform(complex[,] data, bool forward = true)
        {
            var real = data.Real(); unsafe { fixed (double* r = real) DCT(r, real.GetLengths(), forward); }
            var imag = data.Imag(); unsafe { fixed (double* r = imag) DCT(r, imag.GetLengths(), forward); }
            return Mt.Tocomplex(real, imag);
        }

        #endregion

        #region Discrete Wavelet Transform

        unsafe static double inner_loopyarray(double* data, int start, int n, double* filter, int m)
        {
            double a = 0;
            int l = Math.Min(m, n - start);
            for (int k = 0; k < l; k++) a += data[start + k - 0] * filter[k];
            for (int k = l; k < m; k++) a += data[start + k - n] * filter[k];
            return a;
        }
        unsafe static void letaddmul_loopyarray(double* data, int start, int n, double* filter, int m, double scalar)
        {
            int l = Math.Min(m, n - start);
            for (int k = 0; k < l; k++) data[start + k - 0] += filter[k] * scalar;
            for (int k = l; k < m; k++) data[start + k - n] += filter[k] * scalar;
        }
        public unsafe static void DwtLevel(double* data, double* result, int n, double* filter0, double* filter1, int m, bool forward)
        {
            int h = n / 2;
            int o0 = (m / -2) + 2 + n;
            int o1 = (m / -2) + 2 + n;
            if (forward)
            {
                for (int i = 0; i < h; i++)
                {
                    result[i + 0] = inner_loopyarray(data, (i * 2 + o0) % n, n, filter0, m);
                    result[i + h] = inner_loopyarray(data, (i * 2 + o1) % n, n, filter1, m);
                }
            }
            else
            {
                Us.Let(result, 0, n);
                for (int i = 0; i < h; i++) letaddmul_loopyarray(result, (i * 2 + o0) % n, n, filter0, m, data[i + 0]);
                for (int i = 0; i < h; i++) letaddmul_loopyarray(result, (i * 2 + o1) % n, n, filter1, m, data[i + h]);
            }
        }
        // NumericalRecipes3.wtn
        // bug-fixed
        unsafe static void DWT_(double* data, int[] nn, double* f0, double* f1, int m, double* temp0, double* temp1, int level, bool forward)
        {
            int ntot = nn.Product();
            Ex.For(level, l =>
            {
                for (int d = nn.Length, step = 1; --d >= 0; )
                {
                    int n = nn[d];
                    int k = (n >> (forward ? l : level - 1 - l)) & ~1;
                    if (k < m) return;
                    for (int i2 = 0; i2 < ntot; i2 += n * step)
                    {
                        for (int i1 = 0; i1 < step; i1++)
                        {
                            Us.LetStartStep(temp0, data, i1 + i2, step, k);
                            DwtLevel(temp0, temp1, k, f0, f1, m, forward);
                            Us.LetStartStep(data, i1 + i2, step, temp1, k);
                        }
                    }
                    step *= n;
                }
            });
        }
        unsafe static void DWT(double* data, int[] nn, Nm.WaveletType type, int level, bool forward)
        {
            double[] filter0;
            switch (type)
            {
                case WaveletType.Daubechies4: filter0 = WaveletDaubechies4; break;
                case WaveletType.Daubechies6: filter0 = WaveletDaubechies6; break;
                case WaveletType.Daubechies12: filter0 = WaveletDaubechies12; break;
                case WaveletType.Daubechies20: filter0 = WaveletDaubechies20; break;
                default: filter0 = ThrowException.Argument<double[]>("n not yet implemented in Daubs"); break;
            }
            int m = filter0.Length;
            var filter1 = New.Array(m, i => (1 - 2 * (i & 1)) * filter0[m - 1 - i]);
            int nmax = nn.Max();
            if (level == -1) level = int.MaxValue;
            Mt.LetMin(ref level, Mt.Log2Floor(nmax / m));
            var temp0 = new double[nmax];
            var temp1 = new double[nmax];
            fixed (double* f0 = filter0, f1 = filter1, t0 = temp0, t1 = temp1)
            {
                DWT_(data, nn, f0, f1, m, t0, t1, level, forward);
            }
        }
        static readonly double[] WaveletDaubechies4 = new double[4]{
            	+0.4829629131445341, +0.8365163037378079, +0.2241438680420134, -0.1294095225512604
            };
        static double[] WaveletDaubechies6 = new double[] {
                +0.332670552950082616, +0.806891509311092576, +0.459877502118491570, -0.135011020010254589, -0.0854412738820266617, +0.0352262918857095366
            };
        static readonly double[] WaveletDaubechies12 = new double[12]{
            	+0.111540743350, +0.494623890398, +0.751133908021, +0.315250351709, -0.226264693965, -0.129766867567,
            	+0.097501605587, +0.027522865530, -0.031582039318, +0.000553842201, +0.004777257511, -0.001077301085
            };
        static readonly double[] WaveletDaubechies20 = new double[20]{
            	+0.026670057901, +0.188176800078, +0.527201188932, +0.688459039454, +0.281172343661, -0.249846424327,
            	-0.195946274377, +0.127369340336, +0.093057364604, -0.071394147166, -0.029457536822, +0.033212674059,
	            +0.003606553567, -0.010733175483, +0.001395351747, +0.001992405295, -0.000685856695, -0.000116466855,
	            +0.000093588670, -0.000013264203
            };

        public enum WaveletType
        {
            Daubechies4, Daubechies6, Daubechies12, Daubechies20
        }
        public static double[] InverseDiscreteWaveletTransform(double[] data, Nm.WaveletType type, int level = -1) { return DiscreteWaveletTransform(data, type, level, false); }
        public static double[] DiscreteWaveletTransform(double[] data, Nm.WaveletType type, int level = -1, bool forward = true)
        {
            var result = data.CloneX();
            unsafe { fixed (double* r = result) DWT(r, result.GetLengths(), type, level, forward); }
            return result;
        }
        public static double[,] InverseDiscreteWaveletTransform(double[,] data, Nm.WaveletType type, int level = -1) { return DiscreteWaveletTransform(data, type, level, false); }
        public static double[,] DiscreteWaveletTransform(double[,] data, Nm.WaveletType type, int level = -1, bool forward = true)
        {
            var result = data.CloneX();
            unsafe { fixed (double* r = result) DWT(r, result.GetLengths(), type, level, forward); }
            return result;
        }
        public static double[, ,] InverseDiscreteWaveletTransform(double[, ,] data, Nm.WaveletType type, int level = -1) { return DiscreteWaveletTransform(data, type, level, false); }
        public static double[, ,] DiscreteWaveletTransform(double[, ,] data, Nm.WaveletType type, int level = -1, bool forward = true)
        {
            var result = data.CloneX();
            unsafe { fixed (double* r = result) DWT(r, result.GetLengths(), type, level, forward); }
            return result;
        }

        public static complex[] DiscreteWaveletTransform(complex[] data, Nm.WaveletType type, int level = -1, bool forward = true)
        {
            var real = data.Real(); unsafe { fixed (double* r = real) DWT(r, real.GetLengths(), type, level, forward); }
            var imag = data.Imag(); unsafe { fixed (double* r = imag) DWT(r, imag.GetLengths(), type, level, forward); }
            return Mt.Tocomplex(real, imag);
        }
        public static complex[,] DiscreteWaveletTransform(complex[,] data, Nm.WaveletType type, int level = -1, bool forward = true)
        {
            var real = data.Real(); unsafe { fixed (double* r = real) DWT(r, real.GetLengths(), type, level, forward); }
            var imag = data.Imag(); unsafe { fixed (double* r = imag) DWT(r, imag.GetLengths(), type, level, forward); }
            return Mt.Tocomplex(real, imag);
        }

        #endregion

        #endregion

        #region Numerical optimization functions

        // NumericalRecipes3.rtbis
        public static double rtbis(Func<double, double> func, double x1, double x2, double xacc)
        {
            const int maxIteration = 50;
            double dx, xmid, rtb;
            double f = func(x1);
            double fmid = func(x2);
            if (f * fmid >= 0) ThrowException.WriteLine("Root must be bracketed for bisection in rtbis");
            if (f < 0) { dx = x2 - x1; rtb = x1; }
            else { dx = x1 - x2; rtb = x2; }
            for (int j = 0; j < maxIteration; j++)
            {
                fmid = func(xmid = rtb + (dx *= 0.5));
                if (fmid <= 0) rtb = xmid;
                if (Math.Abs(dx) < xacc || fmid == 0) return rtb;
            }
            ThrowException.WriteLine("Too many bisections in rtbis");
            return double.NaN;
        }

        // NumericalRecipes3.rtflsp
        public static double rtflsp(Func<double, double> func, double x1, double x2, double xacc)
        {
            const int maxIteration = 30;
            double xl, xh, del;
            double fl = func(x1);
            double fh = func(x2);
            if (fl * fh > 0) ThrowException.WriteLine("Root must be bracketed in rtflsp");
            if (fl < 0) { xl = x1; xh = x2; }
            else { xl = x2; xh = x1; Mt.Swap(ref fl, ref fh); }
            double dx = xh - xl;
            for (int j = 0; j < maxIteration; j++)
            {
                double rtf = xl + dx * fl / (fl - fh);
                double f = func(rtf);
                if (f < 0.0) { del = xl - rtf; xl = rtf; fl = f; }
                else { del = xh - rtf; xh = rtf; fh = f; }
                dx = xh - xl;
                if (Math.Abs(del) < xacc || f == 0) return rtf;
            }
            ThrowException.WriteLine("Maximum number of iterations exceeded in rtflsp");
            return double.NaN;
        }

        // NumericalRecipes3.rtsec
        public static double rtsec(Func<double, double> func, double x1, double x2, double xacc)
        {
            const int maxIteration = 30;
            double xl, rts;
            double fl = func(x1);
            double f = func(x2);
            if (Math.Abs(fl) < Math.Abs(f)) { rts = x1; xl = x2; Mt.Swap(ref fl, ref f); }
            else { xl = x1; rts = x2; }
            for (int j = 0; j < maxIteration; j++)
            {
                double dx = (xl - rts) * f / (f - fl);
                xl = rts;
                fl = f;
                rts += dx;
                f = func(rts);
                if (Math.Abs(dx) < xacc || f == 0) return rts;
            }
            ThrowException.WriteLine("Maximum number of iterations exceeded in rtsec");
            return double.NaN;
        }

        // NumericalRecipes3.zriddr
        public static double zriddr(Func<double, double> func, double x1, double x2, double xacc)
        {
            const int maxIteration = 60;
            double fl = func(x1);
            double fh = func(x2);
            if ((fl > 0 && fh < 0) || (fl < 0 && fh > 0))
            {
                double xl = x1;
                double xh = x2;
                double ans = -9.99e99;
                for (int j = 0; j < maxIteration; j++)
                {
                    double xm = 0.5 * (xl + xh);
                    double fm = func(xm);
                    double s = Math.Sqrt(fm * fm - fl * fh);
                    if (s == 0.0) return ans;
                    double xnew = xm + (xm - xl) * ((fl >= fh ? 1 : -1) * fm / s);
                    if (Math.Abs(xnew - ans) <= xacc) return ans;
                    ans = xnew;
                    double fnew = func(ans);
                    if (fnew == 0.0) return ans;
                    if (samesign(fm, fnew) != fm) { xl = xm; fl = fm; xh = ans; fh = fnew; }
                    else if (samesign(fl, fnew) != fl) { xh = ans; fh = fnew; }
                    else if (samesign(fh, fnew) != fh) { xl = ans; fl = fnew; }
                    else ThrowException.WriteLine("never get here.");
                    if (Math.Abs(xh - xl) <= xacc) return ans;
                }
                ThrowException.WriteLine("zriddr exceed maximum iterations");
            }
            else
            {
                if (fl == 0) return x1;
                if (fh == 0) return x2;
                ThrowException.WriteLine("root must be bracketed in zriddr.");
            }
            return double.NaN;
        }

        // NumericalRecipes3.zbrent
        public static double FindRootBrent(Func<double, double> function, double start, double end, double tolerance)
        {
            const int maxIteration = 100;
            double a = start, b = end, c = end, d = 0, e = 0;
            double fa = function(a), fb = function(b), fc = fb;
            if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0)) { ThrowException.Argument("Root must be bracketed in zbrent"); return double.NaN; }
            for (int i = 0; i < maxIteration; i++)
            {
                if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { c = a; fc = fa; e = d = b - a; }
                if (Math.Abs(fc) < Math.Abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
                double tol = 2 * Mt.DoubleEps * Math.Abs(b) + tolerance / 2;
                double xm = (c - b) / 2;
                if (Math.Abs(xm) <= tol || fb == 0) return b;
                if (Math.Abs(e) >= tol && Math.Abs(fa) > Math.Abs(fb))
                {  // 逆二乗補間
                    double s = fb / fa, p, q;
                    if (a == c)
                    {
                        p = 2 * xm * s;
                        q = 1 - s;
                    }
                    else
                    {
                        double t = fa / fc;
                        double r = fb / fc;
                        p = s * (2 * xm * t * (t - r) - (b - a) * (r - 1));
                        q = (t - 1) * (r - 1) * (s - 1);
                    }
                    if (p > 0) q *= -1; else p *= -1;
                    if (2 * p < Math.Min(3 * xm * q - Math.Abs(tol * q), Math.Abs(e * q))) { e = d; d = p / q; }  // 成功
                    else { d = xm; e = d; }  // 二分法
                }
                else { d = xm; e = d; }  // 二分法
                a = b; fa = fb;
                b += (Math.Abs(d) > tol) ? d : (xm >= 0 ? tol : -tol);
                fb = function(b);
            }
            ThrowException.WriteLine("FindRootBrent: cannot find a root");
            return double.NaN;
        }

        // NumericalRecipes3.rtnewt
        public static double rtnewt(Func<double, double> funcd, Func<double, double> funcddf, double x1, double x2, double xacc)
        {
            const int JMAX = 20;
            double rtn = 0.5 * (x1 + x2);
            for (int j = 0; j < JMAX; j++)
            {
                double f = funcd(rtn);
                double df = funcddf(rtn);
                double dx = f / df;
                rtn -= dx;
                if ((x1 - rtn) * (rtn - x2) < 0)
                    ThrowException.WriteLine("Jumped out of brackets in rtnewt");
                if (Math.Abs(dx) < xacc) return rtn;
            }
            ThrowException.WriteLine("Maximum number of iterations exceeded in rtnewt");
            return double.NaN;
        }

        // NumericalRecipes3.rtsafe
        public static double rtsafe(Func<double, double> funcd, Func<double, double> funcddf, double x1, double x2, double xacc)
        {
            const int maxIteration = 100;
            double xh, xl;
            double fl = funcd(x1);
            double fh = funcd(x2);
            if ((fl > 0 && fh > 0) || (fl < 0 && fh < 0))
                ThrowException.WriteLine("Root must be bracketed in rtsafe");
            if (fl == 0) return x1;
            if (fh == 0) return x2;
            if (fl < 0) { xl = x1; xh = x2; }
            else { xh = x1; xl = x2; }
            double rts = 0.5 * (x1 + x2);
            double dxold = Math.Abs(x2 - x1);
            double dx = dxold;
            double f = funcd(rts);
            double df = funcddf(rts);
            for (int j = 0; j < maxIteration; j++)
            {
                if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0) || (Math.Abs(2 * f) > Math.Abs(dxold * df)))
                {
                    dxold = dx;
                    dx = 0.5 * (xh - xl);
                    rts = xl + dx;
                    if (xl == rts) return rts;
                }
                else
                {
                    dxold = dx;
                    dx = f / df;
                    double temp = rts;
                    rts -= dx;
                    if (temp == rts) return rts;
                }
                if (Math.Abs(dx) < xacc) return rts;
                f = funcd(rts);
                df = funcddf(rts);
                if (f < 0) xl = rts;
                else xh = rts;
            }
            ThrowException.WriteLine("Maximum number of iterations exceeded in rtsafe");
            return double.NaN;
        }



        static double samesign(double a, double b) { return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }
        static double mulsign(double a, double b) { return b >= 0 ? a : -a; }

        // NumericalRecipes3.Bracketmethod
        static V3<double, double, double> CalcBracket(double a, double b, Func<double, double> func)
        {
            const double GOLD = 1.618034;
            const double GLIMIT = 100;
            const double TINY = 1e-20;
            var ax = a;
            var bx = b;
            var fa = func(ax);
            var fb = func(bx);
            if (fb > fa)
            {
                Mt.Swap(ref ax, ref bx);
                Mt.Swap(ref fb, ref fa);
            }
            var cx = bx + GOLD * (bx - ax);
            var fc = func(cx);
            while (fb > fc)
            {
                double r = (bx - ax) * (fb - fc);
                double q = (bx - cx) * (fb - fa);
                double u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * samesign(Math.Max(Math.Abs(q - r), TINY), q - r));
                double ulim = bx + GLIMIT * (cx - bx);
                double fu;
                if ((bx - u) * (u - cx) > 0)
                {
                    fu = func(u);
                    if (fu < fc)
                    {
                        ax = bx; bx = u;
                        fa = fb; fb = fu;
                        break;
                    }
                    else if (fu > fb)
                    {
                        cx = u;
                        fc = fu;
                        break;
                    }
                    u = cx + GOLD * (cx - bx);
                    fu = func(u);
                }
                else if ((cx - u) * (u - ulim) > 0)
                {
                    fu = func(u);
                    if (fu < fc)
                    {
                        bx = cx; cx = u; u += GOLD * (u - cx);
                        fb = fc; fc = fu; fu = func(u);
                    }
                }
                else if ((u - ulim) * (ulim - cx) >= 0)
                {
                    u = ulim;
                    fu = func(u);
                }
                else
                {
                    u = cx + GOLD * (cx - bx);
                    fu = func(u);
                }
                ax = bx; bx = cx; cx = u;
                fa = fb; fb = fc; fc = fu;
            }
            return New.V3(ax, bx, cx);
        }

        // NumericalRecipes3.Golden
        static V2<double, double> MinimizeGolden(Func<double, double> func, double ax, double bx, double cx)
        {
            const double tol = 3e-8;
            const double R = 0.61803399;
            const double C = 1 - R;

            double x1, x2;
            double x0 = ax;
            double x3 = cx;
            if (Math.Abs(cx - bx) > Math.Abs(bx - ax))
            {
                x1 = bx;
                x2 = bx + C * (cx - bx);
            }
            else
            {
                x2 = bx;
                x1 = bx - C * (bx - ax);
            }
            double f1 = func(x1);
            double f2 = func(x2);
            while (Math.Abs(x3 - x0) > tol * (Math.Abs(x1) + Math.Abs(x2)))
            {
                if (f2 < f1)
                {
                    x0 = x1; x1 = x2; x2 = R * x2 + C * x3;
                    f1 = f2; f2 = func(x2);
                }
                else
                {
                    x3 = x2; x2 = x1; x1 = R * x1 + C * x0;
                    f2 = f1; f1 = func(x1);
                }
            }
            if (f1 < f2)
                return New.V2(x1, f1);
            else
                return New.V2(x2, f2);
        }

        // NumericalRecipes3.Brent
        static V2<double, double> MinimizeBrent(Func<double, double> func, V3<double, double, double> bracket)
        {
            //Console.Write("Brent:");
            const double tolerance = 3e-8;
            const int maxIteration = 100;
            const double goldSection = 0.3819660;
            const double ZEPS = Mt.DoubleEps * 1e-3;

            var a = bracket.v0;
            var b = bracket.v2;
            Mt.LetOrder(ref a, ref b);
            double v, w, x; v = w = x = bracket.v1;
            double fv, fw, fx; fv = fw = fx = func(x);
            double d = 0.0, e = 0.0;
            for (int iteration = 0; ; iteration++)
            {
                if (iteration == maxIteration) { ThrowException.WriteLine("MinimizeBrent: too many iterations"); break; }
                var m = 0.5 * (a + b);
                var tol1 = tolerance * Math.Abs(x) + ZEPS;
                var tol2 = 2.0 * tol1;
                if (Math.Abs(x - m) <= (tol2 - 0.5 * (b - a))) break;

                double be = e, bd = d;
                e = (x >= m ? a : b) - x; d = goldSection * e;  //黄金分割
                if (Math.Abs(be) > tol1)
                {
                    var vw = (x - v) * (fx - fw);
                    var wv = (x - w) * (fx - fv);
                    var p = (x - v) * vw - (x - w) * wv;
                    var q = 2.0 * (vw - wv);
                    if (q > 0) p *= -1; else q *= -1;
                    if (Math.Abs(p) < Math.Abs(0.5 * q * be) && p > q * (a - x) && p < q * (b - x))
                    {
                        e = bd; d = p / q;  //放物線補間
                        if (x + d - a < tol2 || b - (x + d) < tol2) d = mulsign(tol1, m - x);
                    }
                }

                var u = x + (Math.Abs(d) >= tol1 ? d : mulsign(tol1, d));
                var fu = func(u);
                if (fu <= fx)
                {
                    if (u >= x) a = x; else b = x;
                    v = w; fv = fw;
                    w = x; fw = fx;
                    x = u; fx = fu;
                }
                else
                {
                    if (u < x) a = u; else b = u;
                    if (fu <= fw || w == x)
                    {
                        v = w; fv = fw;
                        w = u; fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w)
                    {
                        v = u; fv = fu;
                    }
                }
            }
            //Console.WriteLine();
            return New.V2(x, fx);
        }

        // NumericalRecipes3.Dbrent
        static V2<double, double> MinimizeDbrent(Func<double, double> funcd, Func<double, double> gradient, V3<double, double, double> bracket)
        {
            const double tol = 3e-8;
            const int ITMAX = 100;
            const double ZEPS = Mt.DoubleEps * 1e-3;

            var a = bracket.v0;
            var b = bracket.v2;
            Mt.LetOrder(ref a, ref b);
            double v, w, x; x = w = v = bracket.v1;
            double fv, fw, fx; fw = fv = fx = funcd(x);
            double dv, dw, dx; dw = dv = dx = gradient(x);
            double d = 0.0, e = 0.0;
            for (int iterations = 0; ; iterations++)
            {
                if (iterations == ITMAX) { ThrowException.WriteLine("MinimizeDbrent: too many iterations"); break; }
                var m = 0.5 * (a + b);
                var tol1 = tol * Math.Abs(x) + ZEPS;
                var tol2 = 2.0 * tol1;
                if (Math.Abs(x - m) <= (tol2 - 0.5 * (b - a))) break;
                if (Math.Abs(e) > tol1)
                {
                    var d1 = 2 * (b - a);
                    var d2 = d1;
                    if (dw != dx) d1 = (w - x) * dx / (dx - dw);
                    if (dv != dx) d2 = (v - x) * dx / (dx - dv);
                    var u1 = x + d1;
                    var u2 = x + d2;
                    var ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
                    var ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
                    var olde = e;
                    e = d;
                    if (ok1 || ok2)
                    {
                        if (ok1 && ok2)
                            d = (Math.Abs(d1) < Math.Abs(d2) ? d1 : d2);
                        else if (ok1)
                            d = d1;
                        else
                            d = d2;
                        if (Math.Abs(d) <= Math.Abs(0.5 * olde))
                        {
                            if (x + d - a < tol2 || b - (x + d) < tol2) d = samesign(tol1, m - x);
                        }
                        else
                        {
                            e = (dx >= 0 ? a - x : b - x);
                            d = 0.5 * e;
                        }
                    }
                    else
                    {
                        e = (dx >= 0 ? a - x : b - x);
                        d = 0.5 * e;
                    }
                }
                else
                {
                    e = (dx >= 0 ? a - x : b - x);
                    d = 0.5 * e;
                }

                double u;
                double fu;
                if (Math.Abs(d) >= tol1)
                {
                    u = x + d;
                    fu = funcd(u);
                }
                else
                {
                    u = x + mulsign(tol1, d);
                    fu = funcd(u);
                    if (fu > fx) break;
                }
                var du = gradient(u);
                if (fu <= fx)
                {
                    if (u >= x) a = x; else b = x;
                    v = w; fv = fw; dv = dw;
                    w = x; fw = fx; dw = dx;
                    x = u; fx = fu; dx = du;
                }
                else
                {
                    if (u < x) a = u; else b = u;
                    if (fu <= fw || w == x)
                    {
                        v = w; fv = fw; dv = dw;
                        w = u; fw = fu; dw = du;
                    }
                    else if (fu < fv || v == x || v == w)
                    {
                        v = u; fv = fu; dv = du;
                    }
                }
            }
            return New.V2(fx, x);
        }

        // NumericalRecipes3.Linemethod
        static V2<double[], double> ArgminAlongLine(Func<double[], double> function, double[] argument, double[] direction)
        {
            Func<double, double> newfunction = (double x) => function(Mt.AddMul(argument, direction, x));
            var bracket = CalcBracket(0.0, 1.0, newfunction);
            var brent = MinimizeBrent(newfunction, bracket);
            var arg = Mt.AddMul(argument, direction, brent.v0);
            return New.V2(arg, brent.v1);
        }

        // NumericalRecipes3.Frprmn
        public static Tuple<double[], double> ArgminConjugateGradient(Func<double[], double> function, Func<double[], double[]> gradient, double[] initialArgument, int iterationsMax = 200)
        {
            const double EPS = 1e-18;
            const double toleranceFunction = 3e-8;
            const double toleranceGradient = 1e-8;

            int n = initialArgument.Length;
            var arg = initialArgument.CloneX();
            var func = function(arg);
            var grad = gradient(arg).LetNeg();
            var g = grad.CloneX();
            for (int iterations = 0; ; iterations++)
            {
                if (iterations == iterationsMax)
                {
                    ThrowException.WriteLine("ArgminConjugateGradient: too many iterations");
                    break;
                }
                var h = grad.CloneX();
                var backup = func;
                var result = ArgminAlongLine(function, arg, grad);
                arg = result.v0;
                func = result.v1;
                if (2 * Math.Abs(func - backup) / (Math.Abs(func) + Math.Abs(backup) + EPS) <= toleranceFunction) break;
                grad = gradient(arg).LetNeg();
                if (Ex.Range(n).Max(j => Math.Abs(grad[j]) * Math.Max(Math.Abs(arg[j]), 1)) / Math.Max(func, 1) < toleranceGradient) break;

                var gg = g.SqNorm2();
                if (gg == 0) break;
                double gam = (grad.SqNorm2() - Mt.Inner(grad, g)) / gg;
                g = grad.CloneX();
                grad = Mt.AddMul(g, h, gam);
            }
            return Tuple.Create(arg, func);
        }

        // NumericalRecipes3.Amoeba
        public static Tuple<double[], double> ArgminDownhillSimplex(Func<double[], double> function, double[] initialArgment, double delta, double tolerance, int maxFunctionCall = 5000)
        {
            return ArgminDownhillSimplex(function, initialArgment, New.Array(initialArgment.Length, delta), tolerance, maxFunctionCall);
        }
        public static Tuple<double[], double> ArgminDownhillSimplex(Func<double[], double> function, double[] initialArgment, double[] delta, double tolerance, int maxFunctionCall = 5000)
        {
            var points = New.Array(initialArgment.Length + 1, i =>
            {
                var p = initialArgment.CloneX();
                if (i > 0) p[i - 1] += delta[i - 1];
                return p;
            });
            return ArgminDownhillSimplex(function, points, tolerance, maxFunctionCall);
        }
        public static Tuple<double[], double> ArgminDownhillSimplex(Func<double[], double> function, double[][] initialArgs, double tolerance, int maxFunctionCall = 5000)
        {
            const double TINY = 1e-10;
            var args = initialArgs.Select(v => v.CloneX()).ToArray();
            var argSum = args.Sum();
            var values = args.Select(p => function(p)).ToArray();

            // NumericalRecipes3.Amoeba.amotry
            Func<int, double, double> tryfunc = (int iMax, double fac) =>
            {
                double fac1 = (1.0 - fac) / argSum.Length;
                double fac2 = fac1 - fac;
                var argTry = Mt.Mul(argSum, fac1).LetAddMul(args[iMax], -fac2);
                double valueTry = function(argTry);
                if (valueTry < values[iMax])
                {
                    values[iMax] = valueTry;
                    argSum.LetAdd(argTry).LetSub(args[iMax]);
                    args[iMax] = argTry;
                }
                return valueTry;
            };

            int countFunctionCall = 0;
            int iMin;
            while (true)
            {
                iMin = 0;
                int iMax = values[0] > values[1] ? 0 : 1;
                int iMax2 = 1 - iMax;
                for (int i = 0; i < values.Length; i++)
                {
                    if (values[i] <= values[iMin]) iMin = i;
                    if (values[i] > values[iMax]) { iMax2 = iMax; iMax = i; }
                    else if (values[i] > values[iMax2] && i != iMax) iMax2 = i;
                }
                if (2.0 * Math.Abs(values[iMax] - values[iMin]) / (Math.Abs(values[iMax]) + Math.Abs(values[iMin]) + TINY) < tolerance) break;
                if (countFunctionCall >= maxFunctionCall)
                {
                    //ThrowException.WriteLine("ArgminDownhillSimplex: too many iterations");
                    break;
                }

                countFunctionCall += 2;
                double valueTry = tryfunc(iMax, -1);
                if (valueTry <= values[iMin]) { tryfunc(iMax, 2); continue; }
                if (valueTry < values[iMax2]) { countFunctionCall--; continue; }

                double valueSave = values[iMax];
                valueTry = tryfunc(iMax, 0.5);
                if (valueTry >= valueSave)
                {
                    for (int i = 0; i < values.Length; i++)
                    {
                        if (i == iMin) continue;
                        args[i].LetAdd(args[iMin]).LetMul(0.5);
                        values[i] = function(args[i]);
                    }
                    countFunctionCall += values.Length - 1;
                    argSum = args.Sum();
                }
            }
            return Tuple.Create(args[iMin].CloneX(), values[iMin]);
        }

        // NumericalRecipes1.Ameba
        public static double[] ArgminDownhillSimplex1(Func<double[], double> func, double[] initialArgument, double[] delta, double tolerance)
        {
            const int maxFunctionCall = 5000;
            int countFunctionCall = 0;
            int D = initialArgument.Length;

            double[][] args = New.Array(D + 1, j => New.Array(D, i => initialArgument[i] + (i == j ? delta[i] : 0)));
            double[] argSum = New.Array(D, i => Mt.Sum(D + 1, j => args[j][i]));
            double[] argTry = new double[D];
            double[] values = args.Select(arg => func(arg)).ToArray();

            Func<int, double, double> tryfunc = (int i, double fac) =>
            {
                double[] arg = args[i];
                double fac1 = (1 - fac) / D;
                double fac2 = fac1 - fac;
                for (int j = D; --j >= 0; )
                    argTry[j] = argSum[j] * fac1 - arg[j] * fac2;

                double y = func(argTry);
                if (values[i] > y)
                {
                    values[i] = y;
                    for (int j = D; --j >= 0; )
                    {
                        argSum[j] += argTry[j] - arg[j];
                        arg[j] = argTry[j];
                    }
                }
                return y;
            };

            while (true)
            {
                int iMin = 0;
                int iMax = values[0] > values[1] ? 0 : 1;
                int iMax2 = 1 - iMax;
                for (int i = 0; i < D + 1; i++)
                {
                    double v = values[i];
                    if (v <= values[iMin]) iMin = i;
                    if (v > values[iMax]) { iMax2 = iMax; iMax = i; }
                    else if (v > values[iMax2] && i != iMax) iMax2 = i;
                }

                double d = 2.0 * Math.Abs(values[iMax] - values[iMin])
                    / (Math.Abs(values[iMax]) + Math.Abs(values[iMin]) + Mt.DoubleEps);
                if (d < tolerance) break;
                if (countFunctionCall >= maxFunctionCall)
                {
                    ThrowException.WriteLine("ArgminSimplex: too many function estimations");
                    break;
                }

                countFunctionCall += 2;
                double valueTry = tryfunc(iMax, -1);
                if (valueTry <= values[iMin]) { tryfunc(iMax, 2); continue; }
                if (valueTry < values[iMax2]) { countFunctionCall--; continue; }

                double valueSave = values[iMax];
                valueTry = tryfunc(iMax, 0.5);
                if (valueTry < valueSave) continue;

                countFunctionCall += D;
                for (int i = 0; i < D + 1; i++)
                {
                    if (i == iMin) continue;
                    for (int j = 0; j < D; j++)
                        args[i][j] = (args[i][j] + args[iMin][j]) / 2;
                    values[i] = func(args[i]);
                }
                for (int i = D; --i >= 0; )
                    argSum[i] = Mt.Sum(D + 1, j => args[j][i]);
            }
            return New.Array(D, i => Mt.Average(D + 1, j => args[j][i]));
        }
        #endregion

        #region Association Tests
        public static Tuple<BigInteger, BigInteger> FisherExactTestInteger(int[,] ctable)
        {
            int lengthY = ctable.GetLength(0);
            int lengthX = ctable.GetLength(1);
            int[] restY = New.Array(lengthY, y => Mt.Sum(lengthX, x => ctable[y, x]));
            int[] restX = New.Array(lengthX, x => Mt.Sum(lengthY, y => ctable[y, x]));

            BigInteger nCasesMultinomialX = Mt.MultinomialInteger(restX);
            BigInteger nCasesThis = Mt.Product(lengthY, i => Mt.MultinomialInteger(Ex.Range(lengthX).Select(x => ctable[i, x])));
            BigInteger nCasesLess = 0;
            {
                Action<int, int, BigInteger> function = null;
                function = (int yy, int xx, BigInteger nCasesPrev) =>
                {
                    int backupY = restY[yy];
                    int backupX = restX[xx];
                    int min = Math.Max(0, backupY - Mt.Sum(xx, x => restX[x]));
                    int max = Math.Min(backupY, backupX);
                    if (min > max) ThrowException.Argument("FisherExactTestInteger: ctable");
                    for (int k = max; k >= min; k--)
                    {
                        restY[yy] = backupY - k;
                        restX[xx] = backupX - k;
                        BigInteger nCases = nCasesPrev / Mt.FactorialInteger(k);
                        if (xx > 0)
                        {
                            function(yy, xx - 1, nCases);
                            continue;
                        }
                        if (yy > 1)
                        {
                            function(yy - 1, lengthX - 1, nCases * Mt.FactorialInteger(restY[yy - 1]));
                            continue;
                        }
                        nCases *= Mt.MultinomialInteger(restX);
                        if (nCases <= nCasesThis) nCasesLess += nCases;
                    }
                    restY[yy] = backupY;
                    restX[xx] = backupX;
                };
                function(lengthY - 1, lengthX - 1, Mt.FactorialInteger(restY[lengthY - 1]));
            }
            var gcd = Mt.GreatestCommonDivisor(nCasesMultinomialX, nCasesLess);
            return Tuple.Create(nCasesLess / gcd, nCasesMultinomialX / gcd);
        }

        public static double FisherExactTestDouble(int[,] ctable)
        {
            int lengthY = ctable.GetLength(0);
            int lengthX = ctable.GetLength(1);
            int[] restY = New.Array(lengthY, y => Mt.Sum(lengthX, x => ctable[y, x]));
            int[] restX = New.Array(lengthX, x => Mt.Sum(lengthY, y => ctable[y, x]));

            double nCasesTotal = restX.Sum(c => Mt.LogFactorial(c)) + restY.Sum(c => Mt.LogFactorial(c)) - Mt.LogFactorial(restY.Sum());
            double nCasesThis = ctable.ToEnumerable().Sum(c => Mt.LogFactorial(c)) * (1 - 16 * Mt.DoubleEps);
            double nCasesLess = 0;
            {
                Action<int, int, double> function = null;
                function = (int yy, int xx, double nCasesPrev) =>
                {
                    int backupY = restY[yy];
                    int backupX = restX[xx];
                    int min = Math.Max(0, backupY - Mt.Sum(xx, x => restX[x]));
                    int max = Math.Min(backupY, backupX);
                    if (min > max) ThrowException.Argument("FisherExactTestDouble: ctable");
                    for (int k = max; k >= min; k--)
                    {
                        restY[yy] = backupY - k;
                        restX[xx] = backupX - k;
                        double nCases = nCasesPrev + Mt.LogFactorial(k);
                        if (xx > 0)
                        {
                            function(yy, xx - 1, nCases);
                            continue;
                        }
                        if (yy > 1)
                        {
                            function(yy - 1, lengthX - 1, nCases);
                            continue;
                        }
                        nCases += restX.Sum(c => Mt.LogFactorial(c));
                        if (nCases >= nCasesThis) nCasesLess += Math.Exp(nCasesTotal - nCases);
                    }
                    restY[yy] = backupY;
                    restX[xx] = backupX;
                };
                function(lengthY - 1, lengthX - 1, 0);
            }
            return nCasesLess;
        }

        public static double PearsonChiSquareTest(int[,] ctable)
        {
            int lengthY = ctable.GetLength(0);
            int lengthX = ctable.GetLength(1);
            int[] totalY = New.Array(lengthY, y => Mt.Sum(lengthX, x => ctable[y, x]));
            int[] totalX = New.Array(lengthX, x => Mt.Sum(lengthY, y => ctable[y, x]));
            int total = totalY.Sum();

            double sum = 0;
            for (int y = lengthY; --y >= 0; )
                for (int x = lengthX; --x >= 0; )
                    sum += Mt.Sq(ctable[y, x] - totalY[y] * totalX[x] / (double)total) / (totalY[y] * totalX[x]);
            return Mt.ChiSquareDistributionUpper((lengthY - 1) * (lengthX - 1), sum * total);
        }

        public static double YatesChiSquareTest(int[,] ctable)
        {
            return 1;
        }

        #endregion
    }

    // Human-based functions
    public static partial class Hm
    {
        #region Image functions
        public static double[,] GammaEncoding(this double[,] image, double gamma = 1 / 2.2) { return image.SelectTo(v => Math.Pow(v, gamma)); }
        public static double[,] GammaDecoding(this double[,] image, double gamma = 2.2) { return image.SelectTo(v => Math.Pow(v, gamma)); }

        public static double SrgbEncoding(double v)
        {
            return (v <= 0.00304) ? v * 12.92 : Math.Pow(v, 1 / 2.4) * 1.055 - 0.055;
        }
        public static double SrgbDecoding(double v)
        {
            return (v / 12.92 <= 0.00304) ? v / 12.92 : Math.Pow((v + 0.055) / 1.055, 2.4);
        }
        public static double SrgbEncoding255(double v) { return SrgbEncoding(v) * 255; }
        public static double SrgbDecoding255(double v) { return SrgbDecoding(v / 255); }
        public static double[,] SrgbEncoding(this double[,] image) { return image.SelectTo(v => SrgbEncoding(v)); }
        public static double[,] SrgbDecoding(this double[,] image) { return image.SelectTo(v => SrgbDecoding(v)); }
        public static double[,] SrgbEncoding255(this double[,] image) { return image.SelectTo(v => SrgbEncoding255(v)); }
        public static double[,] SrgbDecoding255(this double[,] image) { return image.SelectTo(v => SrgbDecoding255(v)); }

        public unsafe static double DLSIM(double[,] image0, double[,] image1)
        {
            if (image0.Length != image1.Length) ThrowException.SizeMismatch();
            int n = image0.Length;
            fixed (double* p = image0, q = image1)
            {
                double mx = Us.Sum(p, n) / n;
                double my = Us.Sum(q, n) / n;
                double xx = Us.Sum(p, n, u => Mt.Sq(u - mx)) / n;
                double yy = Us.Sum(q, n, v => Mt.Sq(v - my)) / n;
                double xy = Us.Sum(p, q, n, (u, v) => (u - mx) * (v - my)) / n;
                return 1 - (xy / xx) * (xy / yy);
            }
        }

        public unsafe static double DSSIM(double[,] image0, double[,] image1, double range = 0, double param1 = 0.01, double param2 = 0.03)
        {
            if (image0.Length != image1.Length) ThrowException.SizeMismatch();
            int n = image0.Length;
            double c1 = Mt.Sq(param1 * range);
            double c2 = Mt.Sq(param2 * range);
            fixed (double* p = image0, q = image1)
            {
                double mx = Us.Sum(p, n) / n;
                double my = Us.Sum(q, n) / n;
                double xx = Us.Sum(p, n, x => Mt.Sq(x - mx)) / n;
                double yy = Us.Sum(q, n, y => Mt.Sq(y - my)) / n;
                double xy = Us.Sum(p, q, n, (x, y) => (x - mx) * (y - my)) / n;
                return 1 - (2 * mx * my + c1) * (2 * xy + c2) / ((mx * mx + my * my + c1) * (xx + yy + c2));
            }
        }
        #endregion

        #region Sound functions
        static readonly char[] wavefileChunkRiff = { 'R', 'I', 'F', 'F' };
        static readonly char[] wavefileChunkType = { 'W', 'A', 'V', 'E' };
        static readonly char[] wavefileChunkFrmt = { 'f', 'm', 't', ' ' };
        static readonly char[] wavefileChunkData = { 'd', 'a', 't', 'a' };
        public static void WaveFileSave(string path, double[][] data, int samplesPerSecond, int bitsPerSample = 16)
        {
            try
            {
                using (var writer = new BinaryWriter(File.Create(path)))
                {
                    int channels = data.Length;
                    int dataLength = data[0].Length;

                    int padding = 1;                             // File padding
                    int bytesPerSample = (bitsPerSample / 8) * channels;  // Bytes per sample.
                    int averageBytesPerSecond = bytesPerSample * samplesPerSecond;
                    int chunkDataLength = (bitsPerSample / 8) * dataLength * channels;
                    int chunkFrmtLength = 16;                    // Format chunk length.
                    int chunkRiffLength = chunkDataLength + 36;  // File length, minus first 8 bytes of RIFF description.

                    writer.Write(wavefileChunkRiff);             //4 bytes
                    writer.Write(chunkRiffLength);               //4 bytes
                    {
                        writer.Write(wavefileChunkType);         //4 bytes

                        writer.Write(wavefileChunkFrmt);         //4 bytes
                        writer.Write(chunkFrmtLength);           //4 bytes
                        {
                            writer.Write((short)padding);        //2 bytes
                            writer.Write((short)channels);       //2 bytes
                            writer.Write(samplesPerSecond);      //4 bytes
                            writer.Write(averageBytesPerSecond); //4 bytes
                            writer.Write((short)bytesPerSample); //2 bytes
                            writer.Write((short)bitsPerSample);  //2 bytes
                        }
                        writer.Write(wavefileChunkData);         //4 bytes
                        writer.Write(chunkDataLength);           //4 bytes
                        switch (bitsPerSample)
                        {
                            case 8: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) writer.Write((Byte)Mt.MinMax(data[j][i] * 0x80 + 0x80, 0, 0xff)); break;
                            case 16: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) writer.Write((Int16)Mt.MinMax(data[j][i] * 0x8000, -0x8000, 0x7fff)); break;
                            case 32: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) writer.Write((Int32)Mt.MinMax(data[j][i] * 0x80000000L, -0x80000000L, 0x7fffffffL)); break;
                            default: ThrowException.InvalidOperation("WaveFormat.Save: unknown BitsPerSample"); break;
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.Error.WriteLine(e.Message);
            }
        }
        public static Tuple<double[][], int> WaveFileLoad(string path)
        {
            try
            {
                using (var reader = new BinaryReader(File.OpenRead(path)))
                {
                    if (Ex.Compare(reader.ReadChars(4), wavefileChunkRiff) != 0) ThrowException.InvalidOperation("WaveFormat.Load: not RIFF chunk");
                    reader.ReadInt32();  //chunkRiffLength
                    if (Ex.Compare(reader.ReadChars(4), wavefileChunkType) != 0) ThrowException.InvalidOperation("WaveFormat.Load: not WAVE chunk");

                    int bitsPerSample = 0;
                    int samplesPerSecond = 0;
                    int channels = 0;
                    while (true)
                    {
                        char[] chunkName = reader.ReadChars(4);
                        int chunkLength = reader.ReadInt32();
                        if (Ex.Compare(chunkName, wavefileChunkFrmt) == 0)
                        {
                            reader.ReadInt16();  //shPad                 //2 bytes
                            channels = reader.ReadInt16();               //2 bytes
                            samplesPerSecond = reader.ReadInt32();       //4 bytes
                            reader.ReadInt32();  //averageBytesPerSecond //4 bytes
                            reader.ReadInt16();  //shBytesPerSample      //2 bytes
                            bitsPerSample = reader.ReadInt16();          //2 bytes
                            if (chunkLength > 16) reader.ReadBytes(chunkLength - 16);  //unknown data
                        }
                        else if (Ex.Compare(chunkName, wavefileChunkData) == 0)
                        {
                            int dataLength = chunkLength / channels / (bitsPerSample / 8);
                            var data = New.Array(channels, i => new double[dataLength]);
                            switch (bitsPerSample)
                            {
                                case 8: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) data[j][i] = (reader.ReadByte() - 0x80) / (double)0x80; break;
                                case 16: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) data[j][i] = reader.ReadInt16() / (double)0x8000; break;
                                case 32: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) data[j][i] = reader.ReadInt32() / (double)0x80000000L; break;
                                default: ThrowException.InvalidOperation("WaveFileLoad: unknown BitsPerSample"); break;
                            }
                            return Tuple.Create(data, samplesPerSecond);
                        }
                        else
                        {
                            reader.ReadBytes(chunkLength);  //unknown chunk
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.Error.WriteLine(e.Message);
            }
            return null;
        }
        #endregion
    }

    // Culture-based functions
    public static partial class Ct
    {
        #region Encoding
        // here standard code is utf-8 & \n
        public const string NewLine = "\n";
        public static readonly string[] NewLineCodes = { "\r\n", "\n", "\r" };
        public static readonly Encoding EncodingUTF8 = Encoding.GetEncoding("utf-8");
        public static readonly Encoding EncodingSJIS = Encoding.GetEncoding("shift_jis");
        public static readonly Encoding EncodingEUC = Encoding.GetEncoding("euc-jp");
        public static readonly Encoding EncodingJIS = Encoding.GetEncoding("iso-2022-jp");
        public static readonly Encoding[] Encodings = new[] { EncodingUTF8, EncodingSJIS, EncodingEUC, EncodingJIS };

        public static Encoding DetectEncoding(string path)
        {
            Func<char, bool> IsInvalidChar = c =>
                (c <= 0x1f && c != '\0' && c != '\t' && c != '\n' && c != '\r') || (c == 0x7f) || (c >= 0x80 && c <= 0x9f) || (c >= '\uE000' && c <= '\uF8FF') || (c == '\uFFFD');
            var buffer = new char[4096];
            foreach (var encoding in Encodings)
            {
                using (var file = new StreamReader(path, encoding))
                {
                    while (true)
                    {
                        if (file.EndOfStream) return encoding;
                        var count = file.Read(buffer, 0, buffer.Length);
                        if (buffer.Take(count).Any(c => IsInvalidChar(c))) break;
                    }
                }
            }
            return null;
        }
        #endregion

        #region Words expression
        public static T ConvertFrom<T>(string str)
        {
            return (T)System.ComponentModel.TypeDescriptor.GetConverter(typeof(T)).ConvertFrom(str);
        }

        public static int[] ExtractNumbers(string str)
        {
            var list = new List<int>();
            string[] items = str.Split(',');
            foreach (string s in items)
            {
                if (!s.Contains('-')) { list.Add(int.Parse(s)); continue; }
                string[] t = s.Split('-');
                for (int i = int.Parse(t[0]), j = int.Parse(t[1]); i <= j; ++i)
                    list.Add(i);
            }
            list.Sort();
            return list.ToArray();
        }
        #endregion
    }
}
