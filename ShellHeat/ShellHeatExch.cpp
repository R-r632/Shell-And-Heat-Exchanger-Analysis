#include<iostream>
#include<conio.h>
#include<iomanip>
#include<math.h>
using namespace std;
class Shell
{
    private:
    float T1[10];
    float T2[10];
    float t1[10];
    float t2[10];
    float lmtd[10];
    float eff[10];
    float ntu;
    float q[10];
    float f;
    float u;
    float a;
    float mc;
    float mh;
    int i;
    float cc;
    float c;
    float h;
    float ch;
    float cmin;
    void Calculate()
    {
        for(i=0;i<10;i++)
        {
            lmtd[i]=(((T1[i]-t2[i])-(T2[i]-t1[i]))/log((T1[i]-t2[i])/(T2[i]-t1[i])));
        }
        for(i=0;i<10;i++)
        {
            q[i]=f*u*a*lmtd[i];
        }
        c=mc*cc;
        h=mh*ch;
        if(c<h)
        {
            cmin=c;
        }
        else
        {
            cmin=h;
        }
        ntu=(u*a)/cmin;
        for(i=0;i<10;i++)
        {
            eff[i]=q[i]/(cmin*(T1[i]-t1[i]));
        }
    }
    public:
    Shell()
    {
        f=0.93;
        u=129.7;
        a=0.72;
        mc=2;
        mh=9;
        cc=3.6;
        ch=4.7;
    }
    void Add_Temp()
    {
        cout<<"Enter 10 Entry Temperature Of The Hot Fluid:-"<<endl;
        for(i=0;i<10;i++)
        {
            cout<<"T1["<<i+1<<"]:";
            cin>>T1[i];
        }
        cout<<"Enter 10 Exit Temperature Of The Hot Fluid:-"<<endl;
        for(i=0;i<10;i++)
        {
            cout<<"T2["<<i+1<<"]:";
            cin>>T2[i];
        }
        cout<<"Enter 10 Entry Temperature Of The Cold Fluid:-"<<endl;
        for(i=0;i<10;i++)
        {
            cout<<"t1["<<i+1<<"]:";
            cin>>t1[i];
        }
        cout<<"Enter 10 Exit Temperature Of The Cold Fluid:-"<<endl;
        for(i=0;i<10;i++)
        {
            cout<<"t2["<<i+1<<"]:";
            cin>>t2[i];
        }
    }
    void Given_data()
    {
        cout<<"Correction Factor Based On Type Of Exchanger: "<<f<<endl;
        cout<<"Overall Heat Transfer Coefficient: "<<u<<"W/m^2K"<<endl;
        cout<<"Area Of Inlet and Outlet Of Heat Exchanger: "<<a<<"m^2"<<endl;
        cout<<"Mass per unit second Of The Cold Fluid: "<<mc<<"kg/s"<<endl;
        cout<<"Mass per unit second Of The Hot Fluid: "<<mh<<"kg/s"<<endl;
        cout<<"Specific Heat Capacity Of The Cold Fluid: "<<cc<<"J/kgK"<<endl;
        cout<<"Specific Heat Capacity Of The Hot Fluid: "<<ch<<"J/kgK"<<endl;     
    }
    void Display()
    {
        Calculate();
        cout<<"Log Mean Temperature Difference For The 10 Different Values Entered:-"<<endl;
        for(i=0;i<10;i++)
        {
            cout<<"Log Mean Temperature Difference["<<i+1<<"]:"<<lmtd[i]<<" C"<<endl;
        }
        cout<<"Heat Exchanged For 10 Different Values:-"<<endl;
        for(i=0;i<10;i++)
        {
            cout<<"Heat Exchanged["<<i+1<<"]:"<<q[i]<<" W"<<endl;
        }
        cout<<"Heat Capacity Of The Cold Fluid: "<<c<<" W/K"<<endl;
        cout<<"Heat Capacity Of The Hot Fluid: "<<h<<" W/K"<<endl;
        cout<<"Minimum Heat capacity: "<<cmin<<" W/K"<<endl;
        cout<<"Number of Transfer Units: "<<ntu<<endl;
        cout<<"Effectiveness For The 10 Different Values Entered:- "<<endl;
        for(i=0;i<10;i++)
        {
            cout<<"Effectiveness["<<i+1<<"]:"<<eff[i]<<endl;
        }
    }
};
int main(void)
{
    Shell sh;
    sh.Add_Temp();
    sh.Given_data();
    sh.Display();
}