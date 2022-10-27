#include <iostream>
#include <initializer_list>
#include <cmath>
#include <vector>
using namespace std;
class Matrix
{
public:
    Matrix(Matrix&& c) noexcept :row(c.row), col(c.col), m(c.m)
    {
        /* 带右值引用的拷贝构造，如果等式右端是一个函数返回值，就触发此拷贝构造 */
        c.m = nullptr;
    }
    Matrix(const Matrix& c)
    {
        row = c.row;
        col = c.col;
        m = new double* [row];
        for (int i = 0; i < row; i++)
        {
            m[i] = new double[col];
            for (int j = 0; j < col; j++)
            {
                m[i][j] = c.m[i][j];
            }
        }
    }
    Matrix(int row, int col) :row(row), col(col) //零初始化
    {
        m = new double* [row];
        for (int i = 0; i < row; i++)
        {
            m[i] = new double[col];
            for (int j = 0; j < col; j++)
            {
                m[i][j] = 0;
            }
        }
    }
    Matrix(int row, int col, const initializer_list<double>& x) :row(row), col(col) //赋值初始化
    {
        m = new double* [row];
        auto p = x.begin();
        for (int i = 0; i < row; i++)
        {
            m[i] = new double[col];
            for (int j = 0; j < col; j++)
            {
                if (p != x.end())
                {
                    m[i][j] = *p;
                    ++p;
                }
                else m[i][j] = 0;

            }
        }
    }
    ~Matrix()
    {
        if (m == nullptr) return;
        for (int i = 0; i < row; i++)
        {
            delete[](m[i]);
        }
        delete[]m;
        //cout << "destroy object" << endl;
    }
    void print_Matrix()
    {
        cout << "-----**-----" << endl;
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                cout << " " << m[i][j] << " ";
            }
            cout << endl;
        }
    }
    Matrix& operator+(const Matrix& an_Matrix)
    {
        if ((this->row == an_Matrix.row) && (this->col == an_Matrix.col))
        {
            for (int i = 0; i < row; i++)
            {

                for (int j = 0; j < col; j++)
                {
                    m[i][j] = m[i][j] + an_Matrix.m[i][j];
                }
            }
        }
        else
        {
            //应抛出异常
        }
        return (*this);
    }
    Matrix& operator-(const Matrix& an_Matrix)
    {
        if ((this->row == an_Matrix.row) && (this->col == an_Matrix.col))
        {
            for (int i = 0; i < row; i++)
            {

                for (int j = 0; j < col; j++)
                {
                    m[i][j] = m[i][j] - an_Matrix.m[i][j];
                }
            }
        }
        else
        {
            //应抛出异常
        }
        return (*this);
    }
    Matrix operator*(const Matrix& an_Matrix)
    {
        //注意：不能返回局部变量的引用，因为变量在函数返回时会被释放
        Matrix result(this->row, an_Matrix.col);
        if (this->col == an_Matrix.row)
        {
            for (int i = 0; i < result.row; i++)
            {
                for (int j = 0; j < result.col; j++)
                {
                    //确定result[i][j]的值
                    for (int k = 0; k < col; k++)
                    {
                        result.m[i][j] += m[i][k] * an_Matrix.m[k][j];
                    }
                }
            }
        }
        else
        {
            //应抛出异常
        }
        return result;
    }
    Matrix& operator=(const Matrix& c)
    {
        /* 重载=为深拷贝 */
        row = c.row;
        col = c.col;
        m = new double* [row];
        for (int i = 0; i < row; i++)
        {
            m[i] = new double[col];
            for (int j = 0; j < col; j++)
            {
                m[i][j] = c.m[i][j];
            }
        }
        return (*this);
    }
    Matrix& operator=(Matrix&& c) noexcept
    {
        if (this == &c) return (*this);
        for (int i = 0; i < row; i++)
        {
            //if (m[i]!=nullptr)
            delete[](m[i]); //释放原先的内存
        }
        row = c.row;
        col = c.col;
        m = c.m;
        c.m = nullptr;
        return (*this);
    }
    double* operator[](int row) const
    {
        return m[row];
    }
    Matrix& junior_transform_swap_Row(int row1, int row2)
    {
        /*   初等变换：交换矩阵的两行   */
        if (row1 < row && row2 < row)
        {
            double* p = m[row1];
            m[row1] = m[row2];
            m[row2] = p;
        }
        return (*this);
    }
    Matrix& junior_transform_mul_Row(int row, double x)
    {
        /*   初等变换：某一行乘以一个数   */
        for (int j = 0; j < col; j++) m[row][j] *= x;
        return (*this);
    }
    Matrix& junior_transform_add_Row(int base_row, int ori_row, double x)
    {
        /*   初等变换：ori_row:= ori_row + base_row *x   */
        for (int j = 0; j < col; j++)
        {
            m[ori_row][j] = m[ori_row][j] + m[base_row][j] * x;
        }
        return (*this);
    }
    Matrix& ref_matrix()
    {
        /* 标准阶梯化矩阵 */
        for (int i = 0; i < row - 1; i++)
        {
            //用第i行消去i+1, i+2, ...row行
            //首先寻找第i列，第i, i+1....row行绝对值最大的值，将那一行作为主元
            //交换两个主元，同时交换解的编号
            int max_index = i;
            for (int index = i + 1; index < row; index++)
            {
                if (abs(m[max_index][i]) < abs(m[index][i])) max_index = index;
            }
            junior_transform_swap_Row(i, max_index);
            //cout << "r" << i << "<->" << "r" << max_index << endl;
            //print_Matrix();
            //执行消去
            for (int index = i + 1; index < row; index++)
            {
                double x = m[index][i] / m[i][i];
                junior_transform_add_Row(i, index, -x);
                //cout << "r" << index << ":=" << "r" << index << "-" << x << "*" << "r" << i << endl;
                //print_Matrix();
            }
        }
        //向上消去为单位阵
        for (int i = row - 1; i > 0; i--)
        {
            //用第i行消去i-1, i-2, ...0行
            //执行消去
            for (int index = i - 1; index >= 0; index--)
            {
                double x = m[index][i] / m[i][i];
                junior_transform_add_Row(i, index, -x);
                //print_Matrix();
            }
        }
        for (int i = row - 1; i >= 0; i--) { junior_transform_mul_Row(i, (1 / m[i][i])); }
        //print_Matrix();
        return (*this);
    }

    int getRow() const { return row; }
    int getCol() const { return col; }
protected:
    int row;
    int col;
    double** m;
};
void swap(int* a, int i, int j) { int s = a[i]; a[i] = a[j]; a[j] = s; }
void ref_Matrix(Matrix& m)
{
    /*   高斯主元阶梯化矩阵   */
    int row = m.getRow();
    int col = m.getCol();
    for (int i = 0; i < row - 1; i++)
    {
        //用第i行消去i+1, i+2, ...row行
        //首先寻找第i列，第i, i+1....row行绝对值最大的值，将那一行作为主元
        //交换两个主元，同时交换解的编号
        int max_index = i;
        for (int index = i + 1; index < row; index++)
        {
            if (abs(m[max_index][i]) < abs(m[index][i])) max_index = index;
        }
        m.junior_transform_swap_Row(i, max_index);
        //cout << "r" << i << "<->" << "r" << max_index << endl;
        //m.print_Matrix();
        //执行消去
        for (int index = i + 1; index < row; index++)
        {
            double x = m[index][i] / m[i][i];
            m.junior_transform_add_Row(i, index, -x);
            //cout << "r" << index << ":=" << "r" << index << "-" << x << "*" << "r" << i << endl;
            //m.print_Matrix();
        }
    }
    //向上消去为单位阵
    for (int i = row - 1; i > 0; i--)
    {
        //用第i行消去i-1, i-2, ...0行
        //执行消去
        for (int index = i - 1; index >= 0; index--)
        {
            double x = m[index][i] / m[i][i];
            m.junior_transform_add_Row(i, index, -x);
            //m.print_Matrix();
        }
    }
    for (int i = row - 1; i >= 0; i--) { m.junior_transform_mul_Row(i, (1 / m[i][i])); }
    //m.print_Matrix();
}
class Matrix_E :public Matrix
{
public:
    Matrix_E(int size) :Matrix(size, size)
    {
        for (int i = 0; i < size; i++) m[i][i] = 1;
    }
};
class Col_Vector :public Matrix //列向量矩阵
{
public:
    Col_Vector(int row, const initializer_list<double>& x) :Matrix(row, 1, x) {} //赋值初始化
    Col_Vector(int row) :Matrix(row, 1) {}
    Col_Vector(const Matrix&& m) :Matrix(m) {}
    int get_Length() const { return row; }
};
Matrix cat_matrix_row(const Matrix& m1, const Matrix& m2)
{
    int row1 = m1.getRow();
    int row = row1 + m2.getRow();
    int col = m1.getCol();
    Matrix result = Matrix(row, col);
    /*  以列对齐，将m2粘贴到m1下面行拼接，要求两个矩阵的列数相同,如果不相同，则返回一个零矩阵 */
    if (col == m2.getCol())
    {
        for (int i = 0; i < row1; i++)
        {
            for (int j = 0; j < col; j++)
            {
                result[i][j] = m1[i][j];
            }
        }
        for (int i = row1; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                result[i][j] = m2[i - row1][j];
            }
        }
    }
    else {/*抛出异常(暂时不会) */ }
    return result;
}
Matrix cat_matrix_col(const Matrix& m1, const Matrix& m2)
{
    /*  以行对齐，将m2粘贴到m1右边列拼接，要求两个矩阵的列数相同,如果不相同，则返回一个零矩阵 */
    int row = m1.getRow();
    int col1 = m1.getCol();
    int col = col1 + m2.getCol();
    Matrix result = Matrix(row, col);
    if (row == m2.getRow())
    {
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                if (j < col1) result[i][j] = m1[i][j];
                else result[i][j] = m2[i][j - col1];
            }
        }
    }
    else {/*抛出异常(暂时不会) */ }
    return result;
}
Matrix slice_matrix(const Matrix& m, int start_row, int end_row, int start_col, int end_col)
{
    /* 从矩阵中摘取子矩阵 */
    Matrix result(end_row - start_row + 1, end_col - start_col + 1);
    for (int i = start_row; i <= end_row; i++)
    {
        for (int j = start_col; j <= end_col; j++)
        {
            result[i - start_row][j - start_col] = m[i][j];
        }
    }
    return result;
}
Matrix inv(const Matrix& a)
{
    Matrix b = cat_matrix_col(a, Matrix_E(a.getRow()));
    b.ref_matrix();
    return slice_matrix(b, 0, b.getRow() - 1, b.getRow(), b.getCol() - 1);

}

double product(const Col_Vector& v1, const Col_Vector& v2)
{
    double result = 0;
    if (v1.getRow() == v2.getRow())
    {
        for (int i = 0; i < v1.getRow(); i++) result += v1[i][0] * v2[i][0];
    }
    else
    {
        //抛出异常
    }
    return result;
}
Matrix solve(const Matrix& a, const Col_Vector& b)
{
    Matrix c = cat_matrix_col(a, b); //子类转换为父类是安全的，父类转换为子类是危险的
    c.ref_matrix();
    c.print_Matrix();
    return slice_matrix(c, 0, c.getRow() - 1, c.getRow(), c.getCol() - 1);
}

Matrix polyfit(const Col_Vector& x, const Col_Vector& y, int n)
{
    n = n + 1;
    int length = x.get_Length();
    Matrix A(n, n);
    Col_Vector b(n);
    if (length == y.get_Length())
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i > j)
                {
                    A[i][j] = A[j][i]; //对称阵
                    continue;
                }
                //先生成{1, x,    x^n}在x向量上的序列值
                //i行j列上，多项式次数为i+j
                for (int index = 0; index < length; index++) A[i][j] += pow(x[index][0], i + j);
            }
            for (int index = 0; index < length; index++) b[i][0] += y[index][0] * pow(x[index][0], i);

        }
        A.print_Matrix();
        b.print_Matrix();
        return solve(A, b);
    }

    else
    {
        //抛出异常
    }
    return A;
}

double Romberg(double left, double right, double (*fun)(double), double delta)
{
    /* 龙贝格积分方法 */

    int m_size = 10;
    int i = 1;
    double T = (right - left) / 2 * ((*fun)(left) + (*fun)(right));
    int space_num = 1; //区间个数
    double h = right - left; //区间长度
    Matrix m(m_size, m_size);//T矩阵
    m[0][0] = T;
    double dT = abs(m[i][0] - m[i - 1][0]);
    while (dT > delta and i <= 10)
    {

        if (i + 1 == m_size) //扩容
        {
            Matrix new_m(m_size * 2, m_size * 2);
            for (int a = 0; a < m_size; a++)
            {
                for (int b = 0; b < m_size; b++)
                {
                    new_m[a][b] = m[a][b];
                }
            }
            m_size *= 2;
            m = new_m;
        }
        double U = 0.0;
        cout << "h=" << h << endl;
        for (int index = 0; index < space_num; index++)
        {
            /* 先计算m[0][i]积分，分为2^m[0][i]个区间，每个区间长度h为
            (right-left)/(2^)m[0][i],取（2^m[0][i]）个中点值 */
            U += (*fun)(left + (index + 0.5) * h);
            cout << "x=" << (index + 0.5) * h << endl;
        }
        U *= h;
        cout << "U=" << U << endl;
        T = 0.5 * (T + U); //递推T(2n)

        m[0][i] = T;
        space_num *= 2; //更新区间个数和区间长度
        h /= 2;
        for (int row = 1; row <= i; row++) //计算右上大小为j+1的Romberg(i)三角形
        {
            /* 迭代，直至得到m[i][0] */
            int col = i - row;
            m[row][col] = (pow(4, row) * m[row - 1][col + 1] - m[row - 1][col]) / (pow(4, row) - 1);
        }
        m.print_Matrix();
        dT = abs(m[i][0] - m[i - 1][0]);
        i += 1;
    }
    m.print_Matrix();
    return m[i][0];
}
double test_fun(double x)
{
    //return x*x;
    //return log10(x);
    //return 1/x;
    return exp(x);
}
int main()
{
    //Matrix a(4, 4, { 1, 2, 3, 4,
    //                5, 6, 7, 8 });
    //a[2][0] = 9;
    //a[3][3] = 1;
    //Matrix b(4, 4, { 1, 4, 2, 5, 3, 6 });
    //a.print_Matrix();
    //b.print_Matrix();
    //(a+b).print_Matrix(); //加法测试
    //(a-b).print_Matrix();
    //(a * b).print_Matrix();
    //a.swap_Row(0, 1);
    //a.print_Matrix();
    Matrix a(4, 4, { 2, 2, 3, 6,
                    4, 1, 7, 1,
                    3, 7, 8, 5,
                    3, 5, 7, 9 });
    //    a.print_Matrix();
    //    //ref_Matrix(a);
    //    //矩阵拼接测试
    //    Matrix b = cat_matrix_col(a, Matrix_E(a.getRow()));
    //    //主元标准化测试
    //    b.ref_matrix();
    //    b.print_Matrix();
    //    //矩阵分割提取测试
    //    Matrix c = slice_matrix(b, 0, b.getRow()-1, 0, b.getRow()-1);
    //    c.print_Matrix();
    //    //矩阵求逆测试
    //    Matrix d = inv(a);
    //    d.print_Matrix();
        //cat_matrix_row(a, Matrix(1, 5)).print_Matrix();
        //cat_matrix_col(a, Matrix(4, 1, {1, 2, 3, 4})).print_Matrix();
    //    Col_Vector v(4, {1, 2, 3, 4});
    //    //Matrix e = cat_matrix_col(a, v);
    //    //e.print_Matrix();
    //    Matrix g  = solve(a, v);
    //    g.print_Matrix();
    //    Col_Vector v1 = v;
    //    cout<<product(v, v1)<<endl;

        //v.print_Matrix();
    //    Col_Vector v2(2, {2, 3}), v3(2, {3, 4});
    //    polyfit(v2, v3, 1).print_Matrix();
    Col_Vector v4(7, { 3, 4, 5, 6, 7, 8, 9 }), v5(7, { 2.01, 2.98, 3.50, 5.02, 5.47, 6.02, 7.05 });
    polyfit(v4, v5, 1).print_Matrix();

    //double (*p)(double) = test_fun;
    //cout << (*p)(5) << endl;
    //double d = Romberg(1, 5, p, 10e-8);
    //printf("%f", d);
}
