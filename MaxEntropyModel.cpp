/*In this program,i will show you some details about Derivatives and Partial Derivatives
* In my opinion,program is the appear again about the mathmatic.
*/

#include <iostream>
#include <vector>
#define NUMOFPARAMETER 2
using namespace std;

//该函数用于表达式定义以及返回该点的函数值
double expressionValue(double x)
{
    return pow(x, 2) -7 * x + 12;
}


//计算导数
//第一个参数用于给出x变量，也就是对应的那个点的导数
//第二个参数用于给出Δx的取值(也就相当于导数的精确度了)。
//第二个函数的精度我给到了10的-10次，应该是比较小了，所以如果没有特殊需求的话可以不动。
//因为10的-10次已经非常的小了，近似可以看成是-∞小了。
//之所以给的是pair类型的，因为如果只给double类型的返回值，如果else的返回值和函数的导数值正好相等了，那么就会无法分辨
//所以此处我们保险起见，用一个pair类型的来存储，第二个值用bool类型的，作为一个标识符标识是否计算出了导数。
pair<double,bool> calDerivative(double x,double precision = 1e-10)
{

    //导数存在的充分必要条件：左导数存在且右导数存在而且两者相等。
    //这里我们假设所给的函数是连续的，如果不是连续则要通过手动计算，该方法不可使用
    //可导一定连续，连续不一定可导


    double fValue = expressionValue(x);//计算该点处函数值
    double rightDerivative = (expressionValue(x + precision) - fValue) / precision;//计算右导数
    double leftDerivative = (expressionValue(x - precision) - fValue) / -precision;//计算左导数


    //如果左导等于右导，那么函数可导，控制一下左导和右导的差值范围，因为计算机无法真正的逼近那个值
    //所以我们稍微控制一下误差的范围大小即可，在该范围内，就认定为相等，否则，视为不等。
    //就我看来，10的-6次已经足够小了。所以这里我取10的-6次方，左导数和右导数都有误差，取2倍即可。
    if (abs(rightDerivative - leftDerivative) < 2 * 1e-6)
    {
        return pair<double,bool>(leftDerivative,true);//返回左右导数之中的任意一个都可以
    }
    else
    {
        return pair<double, bool>(leftDerivative, false);//连续但不可导，标识符返回false，表示该点导数不存在。
    }
}

/*
牛顿迭代法(泰勒一阶展开表述)
牛顿迭代法的原理在于用切线方程去求方程的根
具体的原理如下：
我们过函数f(x)的某一点作一条切线，那么该切线方程为：y - y0 = f'(x)(x - x0)
有了上面了导数求法，我们也有y0(也即f(x))，而我们要求的是与x轴的交点,那么此时的y = 0
于是我们就可以通过切线方程直接确定这条切线和x轴的交点，然后再求出这个交点对应的f(x)的值
看看这个f(x)和0的误差，如果误差是在我们接受的范围内，就返回这个x
否则将这个x对应的函数值f(x)所在的那条切线再求出来，按照上面说的步骤继续搞下去(即回到第一步)。
最后根据这个表达式，我们可以得到下一个x的点的表示方法：x = x0 - f(x)/f'(x)
此即为牛顿迭代法，现在我们用牛顿迭代法来求这个方程的根。
此处采用pair<>类型的原因和上面的类似
但是有个弊端，一次只能逼近一个。
*/
pair<double,bool> NewTon_IterateFunc(double x,double precision = 1e-10)//第一个参数随便给就行，它会自动迭代
{
    double fValue = expressionValue(x);//计算x处的函数值
    pair<double,bool> derivative = calDerivative(x);//计算x处的导数值
    if (derivative.second == false) {//异常处理
        return pair<double, bool>(-1, false);
        cout << "导数的求解有问题，请检查精度范围或者是该函数不可导" << endl;
    }


    double newX = x - fValue / derivative.first;//表达式
    int tag = 0;//计算迭代次数，这里我们计算最多50次迭代的结果，如果还是没有，就看看到底是算法的问题还是真的没有解
    while (abs(expressionValue(newX)) > precision)//满足精度条件，就迭代，否则就跳出返回，注意是abs()，要小心负数的条件
    {
        if (tag >= 20) return pair<double, bool>(newX, false);
        newX = newX - fValue / derivative.first;//表达式,注意不要写成了x - fValue / derivative.first，必须是上一次的
        fValue = expressionValue(newX);//计算新的f(x)的值
        derivative = calDerivative(newX);//计算新的x处的导数值
        tag++;
    }

    return pair<double, bool>(newX,true);//返回这个根
}


/*上面是针对一元函数的求导法则，但是对于多元函数的求导法则就是多元函数求偏导了，这时候我们用一个矩阵来表示
这个矩阵我们称为海塞矩阵(我看我们的高等代数的书上好像叫黑塞矩阵)。
* 当目标函数是二次函数时，海塞矩阵退化成一个常数矩阵，从任一初始点出发，牛顿法可一步到达，因此它是一种具有二次收敛性的算法。
对于非二次函数，若函数的二次性态较强，或迭代点已进入极小点的邻域，则其收敛速度也是很快的，这是牛顿法的主要优点。
但是牛顿法的迭代公式中由于没有步长因子，是定步长迭代，对于非二次型目标函数，有时会使函数值上升，即出现f(Xk + 1) > f(Xk)的情况
这是因为没有设置步长，有可能会产生直接越过那个点的情况。
更甚者，可能出现迭代点列{Xk}发散而导致计算失败的情况。
为解决这个问题，出现了“阻尼牛顿法”，增加一个步长因子⋋k，将算法流程的计算公式修改为：
Xk+ 1 = Xk - ⋋k * Hk-1(即黑塞矩阵的逆矩阵) * gk(偏导数在X的第k次迭代后在Xk处所对应的值)
因此，我们需要偏导数的求导描述方法以及矩阵求逆，这个过程还可以优化，即可以不求逆矩阵，因为逆矩阵的计算量实在太大
但是，这里我们还是描述一下这整个过程，即拟牛顿法的迭代过程。
首先，先是对多元函数求偏导的一个过程，比如给x求偏导：
fx(x,y) = (f(x0 + Δx,y0) - f(x0,y0)) / Δx,类似的y的表达式:
fx(x,y) = (f(x0,y0 + Δy) - f(x0,y0)) / Δy，以此类推......
我们描述一下这个过程：
*/


//第一步：首先是给出表达式，表达式必须要自己录入(每个函数的表达式都不一样，实在不好描述)
//参数我们在函数参数中给予填写，这里我以二元函数为例
//同时提供一个重载方法，用来对某一个参数求偏导
//参数为要求偏导的参数是第几个参数
//因为要用到偏导数的表达式，如果不想调参数，就用结构体
//这里我是用结构体(类)来描述。


typedef class MultiFuncDerivative
{
private:
    //vector<pair<double,double>> allParameters

    pair<double,double> x;//输入变量x参数前的常数和x本身的值
    pair<double,double> y;//输入变量y参数前的常数和y本身的值

    //输出表达式z = ?(z代表多元函数的表达式)
    //其pair的第一个参数代表的是x参数前面的常数，第二个参数代表的是x本身是多少
    //用一个容器来描述所有的变量，这样就可以不用实时修改参数的时候还要修改偏导数是对谁求的(也就是写多个if判断)
    vector<pair<double,double>> z;
public:
    //补充一点，如果后面参数过多不好打也可以再次修改程序，通过一个宏变量定义参数个数
    //然后让pair<double,double>进行那个次数的循环，并把数值放入了对应的z里面
    //之后直接把这个步骤封装成一个方法就可以了。
    //那么此时就不需要上面的x和y变量了，取而代之是一个参数数组：
    //vector<pair<double,double>> allParameters
    //哎，既然都写到这了，兄弟们应该也懂我意思了，如果还不了解，可以再来问我。
    //这里我暂时不这样写，因为保持程序的易读性才是我的目的。
    MultiFuncDerivative(pair<double,double> x, pair<double,double> y)
    {
        this->x = x;
        this->y = y;
        z.push_back(x);
        z.push_back(y);
    }


    //计算出Z在对应的x和y处的偏导数值,但是这个没有通用性，我们来个通用性的方法
    //那一个方法只用一个参数就可以了，在这个方法下面，这里我们不采用这个方法。
    //double calZValue(double x, double y)
    //{
    //    z[0].second = x;//把z值对应的点x和y给到z。
    //    z[1].second = y;
    //    zValue += z[0].first * x;
    //    zValue += z[1].first * y;
    //    return zValue;
    //}


    //重载：
    //第一个参数代表的是第几个参数，第二个参数代表的是那个参数是多少
    double calZValue(int i,double parameterVal)
    {
        return z[i].first * parameterVal;
    }
    

    //计算出Z在对应的x和y处的偏导数值,这里x，y传进来的形式我们用vetcor容器，方便操作
    pair<double,bool> expressionOfMultiFunc(vector<double> parametersVal,int num,double precision = 1e-10)
    {
        double zCurVal = 0.0;//计算当前z的值
        double zLeftValueIncremental = 0.0;//计算Z的左增量
        double zRightValueIncremental = 0.0;//计算z的右增量
        //过程：
        for (int i = 0; i < z.size(); ++i)
        {
            double temp = calZValue(i, parametersVal[i]);
            if (num != i)
            {
                zLeftValueIncremental += temp;
                zRightValueIncremental += temp;
                zCurVal += temp;
            }
            else
            {
                zRightValueIncremental += calZValue(i, parametersVal[i] + precision);
                zLeftValueIncremental += calZValue(i, parametersVal[i] - precision);
                zCurVal += temp;
            }
        }


        //计算对应位置的偏导数
        double leftPartialDerivative = (zLeftValueIncremental - zCurVal) / -precision;//左偏导
        double rightPartialDerivative = (zRightValueIncremental - zCurVal) / precision;//右偏导
        if (abs(leftPartialDerivative - rightPartialDerivative) < 2 * precision)
        {
            return pair<double,bool>(leftPartialDerivative,true);
        }
        else
        {
            return pair<double, bool>(-1, false);
        }


        //这里虽然写完了，但其实是非常有限的，什么意思？
        //也就是说这是针对x，y这两个参数而且没有指数的形式
        //所以你会发现数学表达式你必须要自己录入才合适，否则你每次都要写一个新的程序
        //有一个比较完美的方法，就是用结构体来描述x，y.....
        //但是也有弊端，比如我给一个logx，那你这个方法就又没有用了，所以我的建议是自己去写函数表达式就好了
        //同时如果求二阶偏导数的话，计算量就更加的复杂了
        //所以这也是我为什么说最好是自己描述函数就好了
        //这里我写这个偏导数的过程其实就是进行一次描述过程，如果之后函数表达式定下来了
        //那么就可以直接套用这个模板去搞就好了。
    }
}MFD;



int main()
{
    //pair<double, bool> test = NewTon_IterateFunc(6);
    //if (test.second != false) cout << test.first << endl;
    //return 0;

    //Test yes
    //pair<double, double> x = { 2,3 };
    //pair<double, double> y = { 3,5 };
    //MFD* temp = new MFD(x,y);
    //vector<double> a = { 3,5 };
    //cout << temp->expressionOfMultiFunc(a, 0).first << endl;
}
