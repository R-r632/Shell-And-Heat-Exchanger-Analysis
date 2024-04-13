import java.util.*;
public class Shell
{
    private double T1;
    private double T2;
    private double t1;
    private double t2;
    private double lmtd;
    private double eff;
    private double ntu;
    private double q;
    private double f;
    private double u;
    private double a;
    private double mc;
    private double mh;
    private int i;
    private double cc;
    private double c;
    private double h;
    private double ch;
    private double cmin;
    public Shell()
    {
        f=0.93;
        u=129.7;
        a=0.72;
        mc=2;
        mh=9;
        cc=3.6;
        ch=4.7;
    }
    public void Add_Data()
    {
        Scanner in=new Scanner(System.in);
        System.out.println("Enter Entry Temperature Of The Hot Fluid:-");
        T1=in.nextDouble();
        System.out.println("Enter Exit Temperature Of The Hot Fluid:-");
        T2=in.nextDouble();
        System.out.println("Enter Entry Temperature Of The Cold Fluid:-");
        t1=in.nextDouble();
        System.out.println("Enter Exit Temperature Of The Cold Fluid:-");
        t2=in.nextDouble();
    }
    public void Given_data()
    {
        System.out.println("Correction Factor Based On Type Of Exchanger: "+f);
        System.out.println("Overall Heat Transfer Coefficient: "+u+"W/m^2K");
        System.out.println("Area Of Inlet and Outlet Of Heat Exchanger: "+a+"m^2");
        System.out.println("Mass per unit second Of The Cold Fluid: "+mc+"kg/s");
        System.out.println("Mass per unit second Of The Hot Fluid: "+mh+"kg/s");
        System.out.println("Specific Heat Capacity Of The Cold Fluid: "+cc+"J/kgK");
        System.out.println("Specific Heat Capacity Of The Hot Fluid: "+ch+"J/kgK");     
    }
    public void Calculate()
    {
        double z=(Math.log((T1-t2)/(T2-t1)));
        lmtd=(((T1-t2)-(T2-t1))/z);
        q=f*u*a*lmtd;
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
        eff=q/(cmin*(T1-t1));
    }
    public void Display()
    {
        System.out.println("Log Mean Temperature Difference["+i+1+"]:"+lmtd+" C");
        System.out.println("Heat Exchanged["+i+1+"]:"+q+" W");
        System.out.println("Heat Capacity Of The Cold Fluid: "+c+" W/K");
        System.out.println("Heat Capacity Of The Hot Fluid: "+h+" W/K");
        System.out.println("Minimum Heat capacity: "+cmin+" W/K");
        System.out.println("Number of Transfer Units: "+ntu);
        System.out.println("Effectiveness["+i+1+"]:"+eff);
    }
    public static void main(String[] args)
    {
        Shell ob=new Shell();
        ob.Add_Data();
        ob.Given_data();
        ob.Calculate();
        ob.Display();
    }
}