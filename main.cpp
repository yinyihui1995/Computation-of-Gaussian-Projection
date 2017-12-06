//
//  main.cpp
//  ComputationOfGaussianProjection
//
//  Created by 尹羿晖 on 2017/12/3.
//  Copyright © 2017年 尹羿晖. All rights reserved.
//

#include <iostream>
#include <math.h>

void PositiveCalculationOfGaussianProjection(double B,double L,double a,double b,double e1,double e2,double x,double y);
void NegativeComputationOfGaussianProjection(double x,double y,double a,double b,double e1,double e2,double B,double L);
void ConversionOfNeighberStripeCoordinates();
void SetEllipsoidParametric(double &a,double &b);
double SetProjectZone(double L,int &nZoneNumber);
double AngleToDegree(double dAngle);
double DegreeToRadian(double dDegree);
double RadianToDegree(double dRadian);
double RadianToAngle(double dRadian);
void PrintDegree(double dDegree);



int main(int argc, const char * argv[]) {
    
    double B = 0.0,L = 0.0;     //大地坐标 纬度B 经度L
    double x = 0.0,y = 0.0;     //平面坐标 纵坐标x 横坐标y
    double a = 0.0,b = 0.0;   //椭球参数 长半径a 短半径b
    double N;   //卯酉圈曲率半径
    double X;   //当前点到赤道的子午线长度
    double a0,a2,a4,a6,a8;  //计算X的参数
    double m0,m2,m4,m6,m8;  //计算a0,a2,a4,a6,a8的参数
    
    SetEllipsoidParametric(a, b);
    
    //计算椭球的第一、第二偏心率 e1、e2
    double dFirstEccentricityOfEllipsoid = sqrt(pow(a, 2) - pow(b, 2)) / a;
    double dSecondEccentricityOfEllipsoid = sqrt(pow(a, 2) - pow(b, 2)) / b;
    //极点处的子午线曲率半径
    //double c = pow(a, 2) / b;
    
    //选择运算的方式 正算/反算/换带
    std::cout << "选择运算方式 ：1.高斯投影正算 2.高斯投影反算 3.高斯投影邻带换算（未完成）" << std::endl;
    int nOperationMode;
    std::cin >> nOperationMode;
    
    switch (nOperationMode) {
        case 1:
            std::cout << "B" << " " << "L" << std::endl;
            //std::cin >> B >> L;
            B = 30.30;L = 114.20;
            PositiveCalculationOfGaussianProjection(B, L, a, b, dFirstEccentricityOfEllipsoid, dSecondEccentricityOfEllipsoid, x, y);
            //GaussDlg(B, L, x, y);
            break;
        case 2:
            std::cout << "x" << "" << "y" << std::endl;
            //std::cin >> x >> y;
            x = 3380330.773;
            y = 320089.9761;
            NegativeComputationOfGaussianProjection(x, y, a, b, dFirstEccentricityOfEllipsoid,dSecondEccentricityOfEllipsoid, B, L);
            break;
        case 3:
            std::cout << "x" << "" << "y" << std::endl;
            std::cin >> x >> y;
            ConversionOfNeighberStripeCoordinates();
            break;
        default:
            std::cout << "输入错误，请重试" << std::endl;
            
            break;
    }
    
}


//高斯投影正算
void PositiveCalculationOfGaussianProjection(double B,double L,double a,double b,double e1,double e2,double x,double y){
    double M,N;   //子午圈曲率半径 卯酉圈曲率半径
    double X;   //当前点到赤道的子午线长度
    int nZoneNumber; //当前点所在的带号
    double a0,a2,a4,a6,a8;  //用于计算X的参数
    double m0,m2,m4,m6,m8;  //用于计算a0,a2,a4,a6,a8的参数
    double n0,n2,n4,n6,n8;
    //double beta0,beta2,beta4,beta6,beta8;
    
    //double c = pow(a, 2) / b;   //极点子午线曲率半径
    
    double dDegreeB = AngleToDegree(B); //纬度B的十进制形式
    double dRadianB = DegreeToRadian(dDegreeB); //纬度B的弧度形式
    //double dSecondB = dDegreeB * 3600.0;    //以秒为单位的纬度 B”
    
    double dSinB = sin(dRadianB);
    double dCosB = cos(dRadianB);
    //double dCosBSquare = pow(dCosB, 2);
    //double dSinBCosB = sin(dRadianB)*cos(dRadianB);
    
    double t = tan(dRadianB);   //参数 t
    double dEta =sqrt(pow((e2 * dCosB),2));  //参数 𝛈
    double dEtaSquare = pow(e2 * dCosB, 2); // 𝛈^2
    
    double dDegreeL = AngleToDegree(L); //转换成十进制的经度
    int nCentralMeridian = SetProjectZone(L,nZoneNumber);   //中央经线的经度 L0
    nCentralMeridian = 111; //修改中央经线经度为111
    double dDegreeDeltaL = dDegreeL - nCentralMeridian;   //经度L同当前带中央子午线的差值 l
    
    double dDeltaL = DegreeToRadian(dDegreeDeltaL);
    
    //double dSecondDeltaL = dDegreeDeltaL * 3600.0;  //以秒为单位的经差 l“
    
    //double dDeltaLSquare = pow(dDeltaL, 2);
    
    //子午圈曲率半径、卯酉圈曲率半径
    M = a * (1 - pow(e1, 2)) * pow(1 - pow(e1 * dSinB, 2), -3/2);
    N = a / sqrt(1 - pow(e1 * dSinB, 2));
    
    //double dNDeltaLSquare = N * dDeltaLSquare;
    
    m0 = a * (1 - pow(e1, 2));
    m2 = 3 * pow(e1, 2) * m0 / 2;
    m4 = 5 * pow(e1, 2) * m2 / 4;
    m6 = 7 * pow(e1, 2) * m4 / 6;
    m8 = 9 * pow(e1, 2) * m6 / 8;
    
    n0 = a;
    n2 = pow(e1, 2) * n0 / 2;
    n4 = 3 * pow(e1, 2) * n2 / 4;
    n6 = 5 * pow(e1, 2) * n4 / 6;
    n8 = 7 * pow(e1, 2) * n6 / 8;
    
    a0 = m0 + m2 / 2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
    a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7 * m8 / 16;
    a4 = m4 / 8 + 3 * m6 / 16 + 7 * m8 / 32;
    a6 = m6 / 32 + m8 / 16;
    a8 = m8 / 128;
    
    //用带a的式子计算子午线长度
    X = a0 * dRadianB - a2 * sin(2 * dRadianB) / 2 + a4 * sin(4 * dRadianB) / 4 - a6 * sin(6 * dRadianB) / 6 + a8 * sin(8 * dRadianB) / 8;
    
    //X = a0 * dRadianB - sin(dRadianB) * cos(dRadianB) * ((a2 - a4 + a6) + (2 * a4 - 16 * a6 /3) * pow(sin(dRadianB) , 2) + 16 * a6 * pow(sin(dRadianB) , 4) /3);
    
    /*
     //用带beta的式子计算子午线长度
     beta0 = 1 - 3 * pow(e2, 2) / 4 + 45 * pow(e2, 4) / 64 - 175 * pow(e2, 6) / 256 + 11025 * pow(e2, 8) / 16384;
     beta2 = beta0 - 1;
     beta4 = 15 * pow(e2, 4) / 32 - 175 * pow(e2, 6) / 384 + 3675 * pow(e2, 8) / 8192;
     beta6 = - 35 * pow(e2, 6) / 96 + 735 * pow(e2, 8) / 2084;
     beta8 = 315 * pow(e2, 8) / 1024;
     
     X = c * (beta0 * dRadianB +(beta2 * dCosB + beta4 * pow(dCosB, 3) + beta6 * pow(dCosB, 5) + beta8 * pow(dCosB, 7)) * dSinB);
     */
    
    /*
     //用带A的式子计算子午线长度
     double A0 = 1+3.0/4*pow(e1,2)+45.0/64*pow(e1,4)+350.0/512*pow(e1,6)+11025.0/16384*pow(e1,8);
     double A2 = -1.0/2*(3.0/4*e1*e1+60.0/64*pow(e1,4)+525.0/512*pow(e1,6)+17640.0/16384*pow(e1,8));
     double A4 = 1.0/4*(15.0/64*pow(e1,4)+210.0/512*pow(e1,6)+8820.0/16384*pow(e1,8));
     double A6 = -1.0/6*(35.0/512*pow(e1,6)+2520.0/16384*pow(e1,8));
     double A8 = 1.0/8*(315.0/16384*pow(e1,8));
     
     X = a*(1-pow(e1,2))*(A0*dRadianB+A2*sin(2*dRadianB)+A4*sin(4*dRadianB)+A6*sin(6*dRadianB)+A8*sin(8*dRadianB));
     */
    
    //double ax = a0*dRadianB/dDegreeB; //十进制B的系数 a0
    
    //计算x，y
    x = X + N * t * pow(dCosB * dDeltaL, 2) / 2 + N * t * (5 - pow(t, 2) + 9 * dEtaSquare + 4 * pow(dEta, 4)) * pow(dCosB * dDeltaL, 4) / 24 + N * t * (61 - 58 * pow(t, 2) + pow(t, 4)) * pow(dCosB * dDeltaL, 6) / 720;
    
    y = N * dCosB * dDeltaL + N * (1 - pow(t, 2) + dEtaSquare) * pow(dCosB * dDeltaL, 3) / 6 + N * (5 - 18 * pow(t, 2) + pow(t, 4) + 14 * dEtaSquare - 58 * pow(dEta * t, 2)) * pow(dCosB * dDeltaL, 5) / 120;
    
    //double m = cos(dRadianB) * dDeltaL * M_PI / 180.0;
    //x = X + N * t * ( (1/2 + (1/24 * (5 - pow(t,2) + 9 * pow(eta,2) + 4 * pow(eta,4))+ 1 / 720 * (61 - 58 * pow(t,2) + pow(t,4))* pow(m,2))* pow(m,2))*pow(m,2));
    
    //y = N * ((1 + (1/6 * (1 - pow(t,2) + pow(eta,2)) + 1 / 120 *(5 - 18 * pow(t,2) + pow(t,4) + 14 * pow(eta,2) - 58 * pow(eta * t,2))* pow(m,2)) * pow(m,2)) * m);
    
    std::cout << "x:" << x << std::endl;
    std::cout << "y:" << y << std::endl;
    
}


//高斯投影反算
void NegativeComputationOfGaussianProjection(double x,double y,double a,double b,double e1,double e2,double B,double L){
    
    double a0,a2,a4,a6,a8;  //用于计算X的参数
    double m0,m2,m4,m6,m8;  //用于计算a0,a2,a4,a6,a8的参数
    
    m0 = a * (1 - pow(e1, 2));
    m2 = 3 * pow(e1, 2) * m0 / 2;
    m4 = 5 * pow(e1, 2) * m2 / 4;
    m6 = 7 * pow(e1, 2) * m4 / 6;
    m8 = 9 * pow(e1, 2) * m6 / 8;
    
    a0 = m0 + m2 /2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
    a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7 * m8 / 16;
    a4 = m4 / 8 + 3 * m6 / 16 + 7 * m8 / 32;
    a6 = m6 / 32 + m8 / 16;
    a8 = m8 / 128;
    
    //double ax = a0 * M_PI / 180.0;
    
    double X = x;
    
    //迭代法求底点纬度 Bf
    //double Bf = X / a0;   //底点纬度 Bf
    double dRadianBf = X / a0;  //底点纬度 Bf的弧度形式
    
    double Bf_1;
    double FBf;
    double deltaBf;
    //循环迭代部分
    do{
        //FBf = -a2/2 *sin(2*DegreeToRadian(AngleToDegree(Bf))) + a4/4 * sin(4*DegreeToRadian(AngleToDegree(Bf))) - a6/6 * sin(6*DegreeToRadian(AngleToDegree(Bf))) + a8/8 * sin(8*DegreeToRadian(AngleToDegree(Bf)));
        FBf = -a2/2 *sin(2*dRadianBf) + a4/4 * sin(4*dRadianBf) - a6/6 * sin(6*dRadianBf) + a8/8 * sin(8*dRadianBf);
        Bf_1 = (X - FBf)/a0;
        //deltaBf = (RadianToDegree(Bf_1) - RadianToDegree(Bf)) * 3600.0;
        deltaBf =RadianToDegree(Bf_1 - dRadianBf) * 3600.0;
        //Bf = Bf_1;
        dRadianBf = Bf_1;
    }while (fabs(deltaBf) > 1.0e-4);
    
    //double dDegreeBf = AngleToDegree(Bf); //底点纬度 Bf的十进制形式
    //double dRadianBf = DegreeToRadian(dDegreeBf); //底点纬度 Bf的弧度形式
    
    double dDegreeBf = RadianToDegree(dRadianBf);   //底点纬度 Bf的十进制形式
    //double dSecondBf = dDegreeBf * 3600.0;  //以秒为单位的底点纬度 Bf
    
    double tf = tan(dRadianBf);
    
    double etaf = e2 * cos(dRadianBf);  //参数 𝛈f
    double Nf = a / sqrt(1 - pow(e1 , 2) * pow(sin(dRadianBf) , 2));
    
    //double Vf = sqrt(1+pow(e1*cos(dRadianBf), 2));
    double Vf = sqrt(1+pow(etaf , 2));
    
    double dDegreeB = dDegreeBf - 1 / 2 * pow(Vf, 2) * tf * (pow(y/Nf, 2) - 1/12*(5+3*pow(tf, 2)+pow(etaf, 2)-9*pow(etaf*tf, 2))*pow(y/Nf, 4)+1/360*(61+90*pow(tf, 2)+45*pow(tf, 2))*pow(y/Nf, 6))*180/M_PI;
    double dDegreeL = 1/cos(dRadianBf) * (y/Nf - 1/6*(1+2*pow(tf, 2)*pow(etaf, 2))*pow(y/Nf, 3)+1/120*(5+28*pow(tf, 2)+24*pow(tf,2)+6*pow(etaf, 2)+8*pow(etaf*tf, 2))*pow(y/Nf, 5))*180/M_PI;
    
    //std::cout << "B:" << B << std::endl;
    //std::cout << "L:" << L << std::endl;
    
    PrintDegree(dDegreeB);
    PrintDegree(dDegreeL);
    
}


//领带投影坐标换算
void ConversionOfNeighberStripeCoordinates(){
    
    
    
    
    
    
}


//选择椭球体参数
void SetEllipsoidParametric(double &a,double &b){
    std::cout << "选择椭球参数：1.克拉索夫斯基椭球体 2.1975年国际椭球体 3.WGS-84椭球体 4.2000中国大地坐标系（CGCS2000）" << std::endl;
    
    int nEllipsoidParametric;
    std::cin >> nEllipsoidParametric;
    
    switch (nEllipsoidParametric) {
        case 1:
            a = 6378245.0;
            b = 6356863.0187730473;
            break;
        case 2:
            a = 6378140.0;
            b = 6356755.2881575287;
            break;
        case 3:
            a = 6378137.0;
            b = 6356752.3142;
            break;
        case 4:
            a = 6378137.0;
            b = 6356752.3141;
            break;
        default:
            std::cout << "输入错误，请重试" << std::endl;
            SetEllipsoidParametric(a, b);
            break;
    }
}


//选择投影带宽
double SetProjectZone(double L,int &nZoneNumber) {
    int nCentralMeridian; //所在带的中央子午线经度 L0
    int nZoneWide;   //带宽
    double dDegreeL = AngleToDegree(L); //转换成十进制的经度
    std::cout << "选择分带方式：1.3度带 2.6度带" << std::endl;
    std::cin >> nZoneWide;
    //3度带
    if( nZoneWide == 1){
        nZoneNumber = int((dDegreeL + 1.5) / 3);
        nCentralMeridian = 3 * nZoneNumber;
    }
    //6度带
    else {
        nZoneNumber = int(dDegreeL / 6) + 1;
        nCentralMeridian = 6 * nZoneNumber - 3;
    }
    return nCentralMeridian;
}



//将角度（度.分秒）转换成角度（度）
double AngleToDegree(double dAngle){
    int nDegree,nMinute;
    double dSecond,dDegree;
    nDegree = int(dAngle/* + __DBL_EPSILON__*/);
    nMinute = int((dAngle - nDegree) * 100/* + __DBL_EPSILON__*/);
    dSecond = ((dAngle - nDegree) * 100 - nMinute) * 100 ;
    if (dAngle > 0) {
        dDegree = (nDegree + nMinute / 60.0 + dSecond / 3600.0);
    }
    else{
        dDegree = (nDegree - nMinute / 60.0 - dSecond / 3600.0);
    }
    return dDegree;
}

//将角度（度）转换成弧度
double DegreeToRadian(double dDegree){
    return dDegree * M_PI / 180.0;
}

//将弧度转换成角度（度）
double RadianToDegree(double dRadian){
    return dRadian * 180.0 / M_PI;
}

//将弧度转换成角度（度.分秒）
double RadianToAngle(double dRadian){
    int nDegree,nMinute;
    double dSecond,dAngle,dDegree;
    dDegree = RadianToDegree(dRadian);
    nDegree = int(dDegree);
    nMinute = int((dDegree - nDegree) * 60.0);
    dSecond = ((dDegree - nDegree) * 60.0 - nMinute) * 60.0;
    if (dRadian < 0) {
        dAngle = nDegree + nMinute / 100.0 + dSecond / 10000.0;
    }
    else{
        dAngle = nDegree - nMinute / 100.0 - dSecond / 10000.0;
    }
    return dAngle;
}

//将角度（度.分秒）以度分秒（°′″）的形式输出
void PrintDegree(double dDegree){
    int nDegree,nMinute;
    double dSecond;
    nDegree = int(dDegree);
    nMinute = int((dDegree - nDegree) * 60.0);
    dSecond = ((dDegree - nDegree) * 60.0 - nMinute) * 60.0;
    std::cout << nDegree << "°" << nMinute << "′" << dSecond << "″" << std::endl;
}


