//ROS MODEL
// Cortassa 2003 + Generation of free oxygen + ROS scavenging + IMAC (inner membrane anion channel)

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>


using namespace std;

//DEFINE FUNCTIONS

double calculateATPm();
double calculateDmuH();
double calculateNAD();

void sharedVariableUpdate_Mitochondria();

double calculateVCS(); //citrate synthase ACoA,OAA-->CIT
double calculateVACO(); //aconitase CIT->ISOC (isocitrate)
double calculateVIDH(); //isocitrate dehydrogenase ISOC-->aKG
double calculateVKGDH(); //alpha-KG dehydrogenase aKG-->SCoA
double calculateVSL(); //SUC lyase SCoA->SUC

double calculateVSDH(); //SUC-->FUM
double calculateVFH(); //FUM-->MAL
double calculateVMDH(); //MAL-->OAA
double calculateVAAT(); //asparate amino transferase
double calculateVNO();

double calculateVHNe();
double calculateVHFe();

double calculateVhu();
double calculateVATPase();
double calculateVANT();

double calculateVhleak();

double calculateVuni(double Cai);
double calculateVnaca(double Cai);
    

//mito ROS functions
double calculateVIMAC2(); //opening of IMAC
double calculateVt2SO2m();    //transfer of ROS out to cytosol
double calculateVSOD();    //ROS consumption
double calculateVCAT();    //ROS catalase
double calculateVGPX();
double calculateVGR();

//define rate functions
double ATPm, NAD, VCS, VACO, VIDH, VKGDH, VSL, VSDH, VFH, VMDH, VAAT, VNO, VHNe, VHFe, VATPase, Vhu, VANT, Vhleak, VIMAC2, Vt2SO2m, VSOD, VCAT, VGPX, VGR, Vuni, Vnaca;

double dCam, dADPm, dNADH, dISOC, dAKG, dSCoA, dSucc, dFUM, dMAL, dOaa, dROSm, dROSi, dH2O2, dGSH, dDpsi;
double Cam, ADPm, NADH, ISOC, AKG, SCoA, Succ, FUM, MAL, Oaa, ROSm, ROSi, H2O2, GSH, DmuH;
double Dpsi;



//const double Dpsi_low = 0; // Low value
//const double Dpsi_high = 120; // High value
const double Cai_Low = 0.0001; //lower limit is 100nM (0.0001mM)
const double Cai_High = 0.1; //upper limit is 1000nM (0.001mM)
const double pulse_duration = 400; //duration of pulse
const double pulse_frequency = 0.25; //frequency of pulse (Hz)
const double cai_period = 1/pulse_frequency; // Set the period (period is 1/frequency) of the square wave

double AcCoA		= 1.0;          	//1.0AcCoA (mM)
double	b			= 0.5;          	//b
double	Cm			= 1.5;				//Cm (mM) total sum of adenine nucleotides
double	CIK			= 1.0;          	//CIK (mM)
double	Cmito		= 1.812E3;			//Cmito (mM/mV)
                                        //this was 1.812E-3 before

double	CPN			= 10.0;				//CPN (mM)
double	CoA			= 0.02;         	//CoA (mM)
double	DpH			= -0.6;         	//DpH
double	Dpsio		= 50.0;				//Dpsio (mV)* convert from V to mV
double	EtCS		= 0.4;          	//EtCS (mM)

double	EtID		= 0.109;        	//EtID (mM)
double	EtKG		= 0.5;          	//EtKG (mM)
double	EtMD		= 0.154;        	//EtMD (mM)
double	EtSDH		= 0.5;          	//EtSDH (mM)
double	FAD			= 0.01;				//FAD (mM)*** not sure about this

double	FADH2		= 1.24;         	//FADH2(mM)
double	fm			= 0.0003;       	//fm
double	gh			= 1.0E-8;       	//gh (mM/ms/mV)**
double	GLU			= 1.0;         		//GLU (mM)
double	H			= 2.5E-5;       	//H (mM)

double	hm			= 0.5;				//hm
double	KAATeq		= 6.6;          	//KAATeq
double	KaCa		= 0.0005;			//(0.00141)KaCa (mM)
double	KACOeq		= 2.22;         	//KACOeq
double	kact		= 0.00038;			//kact (mM)0.000076	//

double	KADP		= 0.62;				//KADP (mM)
double	Kca			= 0.00127;      	//Kca (mM)
double	kcnsASP		= 1.5E-6; 			//KcnsASP (1/ms)
double	KCS			= 0.5;   			//KCS (1/ms)
double	kfAAT		= 6.44E-4;      	//KfAAT (1/ms)

double	kfACO		= 0.0125;       	//KfACO (1/ms)
double	kfFH		= 0.00332;			//KfFH (1/ms)
double	KFHeq		= 1.0 ;         	//KFHeq
double	kfSL		= 0.005;			//kfSL (1/mM/ms) Is this unit (1/mM/ms) or (1/mM)
double	Kh3			= 6.68E-9;      	//Kh3 (mM)

double	Kh4			= 5.62E-6;      	//Kh4 (mM)
double	Kh1			= 1.131E-5;			//Kh1 (mM)
double	Kh2			= 26.7;         	//Kh2 (mM)
double	kIDH		= 0.05;				//kIDH (1/ms)**
double	KidhNADH	= 0.19;         	//KidhNADH (mM)

double	KiFUM		= 1.3;          	//KiFUM (mM)
double	KiOxaa		= 0.15;         	//KiOxaa (mM)
double	Kioaa		= 0.0031;       	//Kioaa (mM)
double	kKGDH		= 7.5E-2;			//##(7.5e-2)kKGDH (1/ms)(7.50e-2*0.69)
double	KmAcCoA		= 0.0126;			//KmAcCoA (mM)

double	Kmal		= 1.493;        	//Kmal (mM)
double	kMDH		= 0.111;			//KMDH (1/ms)
double	Kmg			= 0.0308;       	//Kmg (mM)
double	KmIDNAD		= 0.923;        	//KmIDNAD (mM)
double	Kmiso		= 1.52;         	//Kmiso (mM)

double	KmKG		= 1.94;         	//KmKG (mM)
double	KmKGNAD		= 38.7 ;        	//KmKGNAD (mM)
double	KmmNAD		= 0.2244;       	//KmmNAD (mM)
double	KmSucc		= 0.03 ;        	//KmSucc (mM)
double	KmOaa		= 0.00064;      	//KmOaa (mM)

double	Koff		= 0.0399;       	//Koff
double	kSDH		= 0.005;			//##(0.005)kSDH (1/ms)(5e-3*1.17)
double	KSLeq		= 3.115;        	//KSLeq
double	Mg			= 0.4;          	//Mg (mM)
double	nID			= 2.0 ;         	//nID

double	nKG			= 1.2;          	//nKG
double	rhoF1		= 0.05; 			//(0.05)rhoF1 (mM)**(0.05*1.27)
double	rhoREN		= 1E-4;				//1.0E-4		//##(1.0e-1)rhoREN (mM)**(1.0e-1*1.17)
double	rhoREF		= 3.75E-4;		//rhoREF (mM)
double	VmDT		= 5e-3;			//VmDT (mM/ms)**

double	VmNC		= 1.0E-4;		//VmNC (mM/ms)**
double	Vmuni		= 0.0275;		//Vmuni (mM/ms)**
double	Pi			= 20.0;			//20 Pi (mM)
double	ATPi 		= 6.25;
double	Nai 		= 10;


//double	Cai			= 0.0001;//make this a square function
double Cai;

double	ADP			= 0.05;

double	af 			= 1.0e4; 		//Activation factor by cytoplasmic O2-
double	Kcc 		= 1.0e-2;		//Activation constant of IMAC by O2-  (mM)
double	G_L 		= 0.035e-6;  	//Leak conductance for IMAC (mM mS-1 mV-1)
double	G_max 		= 3.9085e-6;	//Integral conductance of IMAC at saturation (mM mS-1 mV-1)
double	Kappa 		= 70e-3 ;		//70e-3Steepness factor (mV-1)

double	Em 			= 0.004e3;		//Potential at half saturation (mV)
double	R 			= 8.315e3;		//Gas constant (mV C mol-1 K-1)
double	T 			= 310.16;		//Temperature (K)
double	FF 			= 96480;		//Farady constant (C mol-1)
double	k1_SOD 		= 1.2e3; 		//Second-order rate constant of conversion between native
								// oxidized and reduced SOD and it inactive form (mM-1 mS-1)
double	k5_SOD 		= 0.25e-3; 		//First-order rate constant for conversion between inactive and active oxidized SOD (1/mS)
double	k3_SOD 		= 2.4e1 ;		//2nd-order rate constant of conversion between native reduced SOD and its inactive form (mM/mS)
double  etSOD 		= 1.43e-3;		//1.43e-3 	//Intramitoular concentration of SOD (mM)
double	Kki	 		= 0.5;			//Inhibition constant for H2O2 (mM)
double	k1_CAT 		= 1.7e1; 		//1.7e1Rate constant of CAT (mM-1 mS-1)

double	EtCAT 		= 0.01;	 		//Intracellular concentratio of CAT (mM)
double	fr 			= 0.05;			//Hydrogen peroxide inhibition factor for CAT
double	EtGPX 		= 0.01;	 		//0.0Intracellular conc of GPX (mM)
double	Phi1 		= 0.5e-2; 		//Constant for GPX activity (mM mS)
double	Phi2 		= 0.75e-0;		//Constant for GPX activity (mM mS)

double	k1_GR 		= 5.0e-3; 		//5.0e-7Rate constant of GR (mS-1)
double	EtGR 		= 0.01;	 		//0.01Intracellular conc of GR (mM)
double	Km_GSSG 	= 0.06;	 		//M-M constant for oxidized glutathione of GR (mM)
double	Km_NADPH 	= 0.015 ;		//M-M constant for NADH for GR
double	NADPH 		= 1.0; 			//

double	GT 			= 1.0;			//Total glutathione pool
	//shunt		= 5e-2	;		//18e-2	0.18 shunt
double	j_Vt2SO2m	= 0.1;

double	RT_over_F	= 26.729246592691901;	//RToverF
double	kh_1		= 8.1E-5;		//kh_1 (mM)
double	kh_2		= 5.98E-5;		//kh_2 (mM)
double	ra			= 6.394E-13;	//ra (1/ms)
double	rb			= 1.762E-16;	//rb (1/ms)
double	rc1			= 2.656E-22;	//rc1 (1/ms)

double	rc2			= 8.632E-30;	//rc2 (1/ms)
double	r1			= 2.077E-18;	//r1
double	r2			= 1.728E-9;		//r2
double	r3			= 1.059E-26;	//r3
double	kres		= 1.35E18;		//kres

double	kresf		= 5.765E13;		//kresf
double	g			= 0.85;			//g
double	pa			= 1.656E-8;		//pa (1/ms)
double	pb			= 3.373E-10;	//pb (1/ms)
double	pc1			= 9.651E-17;	//pc1 (1/ms)

double	pc2			= 4.585E-17; 	//pc2 (1/ms)
double	p1			= 1.346E-8;		//p1
double	p2			= 7.739E-7;		//p2
double	p3			= 6.65E-15;		//p3
double	kf1			= 1.71E4;		//kf1

double	ktrans		= 0.019;		//ktrans (mM)0.0038	//
double	L			= 110.0	;		//L
double	na			= 2.8;			//na
double	Kna			= 9.4;			//Kna (mM)
double	Knca		= 3.75E-4;		//Knca, mislabeled: Kca (mM)

double	n			= 3	;			//n


//defining varaibles used in functions (making global so I don't have to redefine below)
double inv_ATPi = 1.0 / ATPi;
double F_over_RT = 1.0 / RT_over_F;
double DmuH_Constant = -2.303 * RT_over_F * DpH;
double VCS_C1 = (KCS * EtCS * AcCoA) / (KmAcCoA + AcCoA);
double one_inv_KACOeq = 1.0 + 1.0 / KACOeq;
    
double inv_KADP = 1.0 / KADP;
double inv_KaCa = 1.0 / KaCa;
double inv_KidhNADH = 1.0 / KidhNADH;
double kIDH_EtID = kIDH * EtID;
double VIDH_Constant = 1.0 + H / kh_1 + kh_2 / H;
    
    
    //Must be done after calculateNAD is Called
    //double KmIDNAD_NAD = KmIDNAD / NAD;
double KmIDNAD_NAD, FRT_6_g, exp_FRT_6_g_DmuH;
    
    //After calculateDmuH
    //double FRT_6_g = 6.0 * g / RT_over_F;
    //double exp_FRT_6_g_DmuH = exp(FRT_6_g * DmuH);
     
    
double Mg_Kmg_1 = Mg / Kmg + 1.0;
double Mg_Kmg_1_Kca = Mg_Kmg_1 / Kca;
double kKGDH_EtKG = kKGDH * EtKG;
double KmKGNAD_KmIDNAD = KmKGNAD / KmIDNAD;
double CoA_KSLeq = CoA / KSLeq;
double kSDH_EtSDH = kSDH * EtSDH;
double KmSucc_KiFUM = KmSucc / KiFUM;
double inv_KiOxaa = 1.0 / KiOxaa;
double kfFH_KFHeq = kfFH / KFHeq;
double kMDH_Fh_EtMD = pow( ( 1.0 / ( 1.0 + Kh3 / H + Kh3 * Kh4 / pow(H,2) ) ) ,2) *
    ( 1.0 / ( 1.0 + H / Kh1 + pow(H,2) / ( Kh1 * Kh2 ) ) + Koff) *
    kMDH * EtMD;
double Kmal_Kioaa = Kmal / Kioaa;
double VAAT_Constant = kfAAT * GLU * kcnsASP * KAATeq / kfAAT;
double kcnsASP_KAATeq_kfAAT = kcnsASP * KAATeq / kfAAT;
double KfAAT_GLU = kfAAT * GLU;
double KfAAT_KAATeq = kfAAT/KAATeq;
double kres_sq_KmIDNAD = kres * kres / KmIDNAD;
double exp_6_FRT_Dpsio = exp( 6.0 * Dpsio / RT_over_F);
double r1_exp_6_FRT_Dpsio = r1 * exp_6_FRT_Dpsio;
double ra_rc1_exp_6_FRT_Dpsio = ra + rc1 * exp_6_FRT_Dpsio;
    
double rhoREN_ra_rc1_exp_6_FRT_Dpsio = 0.5 * rhoREN * ra_rc1_exp_6_FRT_Dpsio;
double rhoREN_rc2 = 0.5 * rhoREN * rc2;
double rhoREN_ra = 0.5 * rhoREN * ra;
double rhoREN_6_ra = 6.0 * rhoREN * ra;
double rhoREN_6_ra_rb = 6.0 * rhoREN * ( ra + rb );
    
double AREF = RT_over_F * log ( kresf * sqrt ( FADH2 / FAD ) );
double exp_AREF_FRT = exp( AREF / RT_over_F );
double ra_exp_AREF_FRT = 4.0 * ra * exp_AREF_FRT;
double ra_rb = 4.0 * (ra + rb);
double VFO_VHFe_C1 = ( 1.0 + r1 * exp_AREF_FRT) * exp_6_FRT_Dpsio;
double r2_r3_exp_AREF_FRT = r2 + r3 * exp_AREF_FRT;
double exp_3_FRT_Dpsio = exp( 3.0 * Dpsio / RT_over_F);
double FRT_3 = 3.0 / RT_over_F;
double kf1_Pi = kf1 / Pi;
double VATPase_C1 = ( 100.0 * pa + pc1 * exp_3_FRT_Dpsio);
double pa_pb_3 = 3.0 * (pa + pb);
double pa_300 = 300.0 * pa;
double p1_exp_3_FRT_Dpsio = p1 * exp_3_FRT_Dpsio;
    
double hm_F_over_RT = hm * F_over_RT;
double VmDT_75 = 0.75 * VmDT;
double VmDT_20 = 20.0 * VmDT;
double tenDiv9 = 10.0 / 9.0;
    
double inv_ktrans = 1.0 / ktrans;
double inv_kact = 1.0 / kact;
double Vmuni_ktrans = Vmuni / ktrans;
double FRT2 = 2.0 * F_over_RT;
double b_05 = b * 0.5;
double FRT2_Dpsi = FRT2 * ( (Dpsi) - 91.0 );
    
double inv_Cmito = 1.0 / Cmito;
double two_b = 2.0 * b;


double I_Total, dV, V;
double I_inj;
double tt, dt, i;

int cycle_length, total_beats, beat;

int main() {

    int counter = 0;
    total_beats = 5;
    cycle_length = 500;
    dV = 0.0;
    dt = 0.005;
    V = -84;
    
    //define output files
    ofstream myfile("outputROS_HZ_DE1.txt");
    ofstream myfile2("outputROS_HZ_V.txt");


    tt = 0.0;
for (tt = 0; tt<= (total_beats*cycle_length); tt=tt+dt) {
    //0.5 is a duration for a stimulus, want to apply a 0.5 duration stimulus every beat

//initial conditions
       Cam = 1.13E-04;
       ADPm = 0.00953;
       Dpsi = 124.498;
       NADH = 2.0803;
       ISOC = 0.10829;
       AKG  = 2.82E-04;
       SCoA = 0.82699;
       Succ = 4.88E-04;
       FUM = 0.01034;
       MAL = 0.00375;
       Oaa = 4.74E-08;
       ROSm = 0.00238;
       ROSi = 1.19E-08;
       H2O2 = 8.82E-09;
       GSH = 0.84055;

    
    double phase = fmod(tt, cai_period);
    //calculate the time within the pulse cycle
    double pulse_time = fmod(phase, cai_period);
    //check if the current time is within the pulse duration
    //if (phase < psi_period /2) { //this is the old line
    if (pulse_time < pulse_duration/1000.0){
        Cai = Cai_High;
    } else {
        Cai = Cai_Low;
    }
     


    ATPm = calculateATPm();
    DmuH = calculateDmuH();
    NAD = calculateNAD();
    
    sharedVariableUpdate_Mitochondria();
    
    VCS = calculateVCS(); //citrate synthase ACoA OAA --> CIT
    VACO = calculateVACO(); //aconitase CIT --> ISOC
    VIDH = calculateVIDH(); //ISOC --> aKG
    VKGDH = calculateVKGDH(); //aKG --> SCoA
    VSL = calculateVSL(); //SCoA --> SUC
    VSDH = calculateVSDH(); //SUC --> FUM
    VFH = calculateVFH(); //FUM --> MAL
    VMDH = calculateVMDH(); //MAL --> OAA
    VAAT = calculateVAAT(); //asparate amino transferase
    VNO = calculateVNO();
    VHNe = calculateVHNe(); //NADH oxidation NADH --> NAD
    VHFe = calculateVHFe(); //FADH2 oxidation FADH2 --> FAD
    Vhu = calculateVhu();
    VATPase = calculateVATPase();
    VANT = calculateVANT();
    Vhleak = calculateVhleak(); //h leakage
    
    //mito ROS
    VIMAC2 = calculateVIMAC2(); //opening of IMAC
    Vt2SO2m = calculateVt2SO2m();	//transfer of ROS out to cytosol
    VSOD = calculateVSOD();	//ROS consumption
    VCAT = calculateVCAT();	//ROS catalase
    VGPX = calculateVGPX();
    VGR = calculateVGR();
    //
    Vuni = calculateVuni(Cai);
    Vnaca = calculateVnaca(Cai);
    
     
    
    //dDpsi = (VHNe+VHFe-Vhu-VANT-Vhleak-Vnaca-2*Vuni)*inv_Cmito;

    dCam  =  fm * (Vuni - Vnaca);
    dADPm =  VANT - VATPase - VSL;
    
    dDpsi = -(-VHNe - VHFe + Vhu + VANT + Vhleak + two_b * Vnaca + 2.0 * Vuni + VIMAC2) * inv_Cmito; //eliminate a few of these at a time to see what is giving NAN
    
    dNADH =  -VNO + VIDH + VKGDH + VMDH;
    dISOC =  VACO - VIDH;
    dAKG  =  VIDH + VAAT - VKGDH;
    dSCoA =  VKGDH - VSL;
    dSucc =  VSL - VSDH;
    dFUM  =  VSDH - VFH;
    dMAL  =  VFH - VMDH;
    dOaa  =  VMDH - VCS - VAAT;
   
    dROSm =  1 * VNO - Vt2SO2m;
    dROSi =  Vt2SO2m - VSOD ;
    dH2O2 =  VSOD - VCAT - VGPX;
    dGSH  =  -VGPX + VGR;


    Cam = Cam + dCam*dt;
    ADPm = ADPm  + dADPm*dt;
    Dpsi = Dpsi + dDpsi*dt;
    NADH = NADH +  dNADH*dt;
    ISOC = ISOC +  dISOC*dt;
    AKG  = AKG +  dAKG*dt;
    SCoA = SCoA +  dSCoA*dt;
    Succ = Succ +  dSucc*dt;
    FUM =FUM +  dFUM*dt;
    MAL = MAL +  dMAL*dt;
    Oaa = Oaa +  dOaa*dt;
    ROSm = Oaa +  dOaa*dt;
    ROSi = ROSi +  dROSi*dt;
    H2O2 = H2O2 +  dH2O2*dt;
    GSH = GSH +  dGSH*dt;
    
    

    counter += 1;

    if (tt>0 && counter%10==0){
        myfile << tt << "," << Dpsi << "," << Cam << "," << ADPm << "," << NADH << "," << ISOC << "," << AKG << "," << SCoA << "," << Succ << "," << FUM << "," << MAL << "," << Oaa << "," << ROSm << "," << ROSi << "," << H2O2 << "," << GSH << "," << ATPm << "," << VSL << endl;
        
        myfile2 << tt << "," << VANT << "," << VATPase << "," << VSL << "," << Cai << "," << Vuni << "," << Vnaca << endl;
    }
}
return 0;

}



//Mitochondia; Intemediate Variable; Optimize Things that use this
double calculateATPm()
{
    ATPm = Cm - (ADPm);
    return ATPm;

}

double calculateDmuH(){
  	DmuH = DmuH_Constant + (Dpsi) ;
return DmuH;

  }


//Mitochondia; Intemediate Variable;  Optimize Things that use this
double calculateNAD()
{
    NAD = CPN - (NADH);
    return NAD;

}

//Mitochondia; Intemediate Variable;  Optimize Things that use this
double calculateVCS()
{
	//Parameter: KmOaa, KCS, EtCS, AcCoA, KmAcCoA
		VCS = VCS_C1 * (Oaa) / ((Oaa) +  KmOaa);
		return VCS;

}

//Mitochondia;Also CIT; Intemediate Variable;  Optimize Things that use this
double calculateVACO()
{
		VACO = kfACO * ( CIK - (AKG) - (SCoA) - (Succ) - (FUM) - (MAL) - (Oaa) - (ISOC)* one_inv_KACOeq );
		return VACO;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this;1/Isoc;Dependent on NAD;KmIDNAD_NAD
double calculateVIDH()
{
		double Fa = 1.0 / (( 1.0 + (ADPm) * inv_KADP) * (1.0 + (Cam) * inv_KaCa));
		double Fi = 1.0 + (NADH) * inv_KidhNADH;
		VIDH = kIDH_EtID /
		  (VIDH_Constant + KmIDNAD_NAD * Fi + pow( Kmiso/(ISOC), nID) * Fa * (1.0 + KmIDNAD_NAD  * Fi));
		return VIDH;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this; Dependent on Calculate NAD(KmIDNAD_NAD); 1/akg
double calculateVKGDH()
{
		double a = ( Mg_Kmg_1 + Mg_Kmg_1_Kca * (Cam) );
		VKGDH = kKGDH_EtKG * a / ( a + pow(KmKG/AKG,nKG) + (KmKGNAD_KmIDNAD * KmIDNAD_NAD) );

		//VKGDH = kKGDH_EtKG * a / ( a + KmKG + (AKG) * pow(KmKGNAD_KmIDNAD * KmIDNAD_NAD,nKG) );
		return VKGDH;
}

//Mitochondia;Also CIT; Intemediate Variable;  Optimize Things that use this; Dependent on calculateATPm
double calculateVSL()
{
		VSL = kfSL * ( (SCoA) * (ADPm) - CoA_KSLeq * (Succ) * ATPm);
		return VSL;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this;1/Succ
double calculateVSDH()
{
		//Parameters: kSDH, EtSDH, KmSucc, KiFUM, KiOxaa
		VSDH = kSDH_EtSDH * (Succ) / ((Succ) + (KmSucc + KmSucc_KiFUM * (FUM)) * (1.0 + inv_KiOxaa * (Oaa)) );
		return VSDH;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this;Dependent on NAD(KmIDNAD_NAD)
double calculateVFH()
{
		//Parameters: kfFH, KFHeq
		VFH = kfFH * (FUM) - kfFH_KFHeq * (MAL);
		return VFH;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this;
double calculateVMDH()
{
		//Parameters: Kh1, Kh2, Kh3, Kh4, Koff, H, Kioaa, KmmNAD, Kmal, EtMD, kMDH
		VMDH = kMDH_Fh_EtMD * (MAL) * NAD /
			( ( (MAL) + Kmal + (Oaa) * Kmal_Kioaa ) * ( KmmNAD + NAD ) );
		return VMDH;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this;
double calculateVAAT()
{
		//Parameters: kfAAT, GLU, kcnsASP, KAATeq,
		VAAT = VAAT_Constant * (Oaa) / ( kcnsASP_KAATeq_kfAAT + (AKG) );
		return VAAT;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this;RT_over_F; F_over_RT;KmIDNAD_NAD;DmuH;exp_6_FRT_Dpsio;exp_FRT_6_g_DmuH;
double calculateVNO()
{		//Parameters: kres, rhoREN, Dpsio, g, ra, rc1, r1, rc2, rb
		double AREN =  sqrt ( (NADH) * kres_sq_KmIDNAD * KmIDNAD_NAD) ;
		double denominator1 = 1.0 / ( (exp_6_FRT_Dpsio + r1_exp_6_FRT_Dpsio *  AREN ) + ( r2 + r3 *  AREN ) * exp_FRT_6_g_DmuH );
		VNO  = ( (rhoREN_ra_rc1_exp_6_FRT_Dpsio + rhoREN_rc2 * exp_FRT_6_g_DmuH) *  AREN  - rhoREN_ra * exp_FRT_6_g_DmuH ) * denominator1;
		return VNO;
		//VHNe = (rhoREN_6_ra * AREN  - rhoREN_6_ra_rb * exp_FRT_6_g_DmuH ) *	denominator1;
}

double calculateVHNe()
{		//Parameters: kres, rhoREN, Dpsio, g, ra, rc1, r1, rc2, rb
		double AREN =  sqrt ( (NADH) * kres_sq_KmIDNAD * KmIDNAD_NAD) ;
		double denominator1 = 1.0 / ( (exp_6_FRT_Dpsio + r1_exp_6_FRT_Dpsio *  AREN ) + ( r2 + r3 *  AREN ) * exp_FRT_6_g_DmuH );

		VHNe = (rhoREN_6_ra * AREN  - rhoREN_6_ra_rb * exp_FRT_6_g_DmuH ) *	denominator1;
		return VHNe;
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;RT_over_F;exp_6_FRT_Dpsio;exp_FRT_6_g_DmuH;
double calculateVHFe()
{		//Parameters: kresf, FADH2, FAD, r2, r3, rhoREF
		double denominator = rhoREF / (VFO_VHFe_C1 + r2_r3_exp_AREF_FRT * exp_FRT_6_g_DmuH);
		VHFe = ( ra_exp_AREF_FRT - ra_rb*exp_FRT_6_g_DmuH ) * denominator;
		return VHFe;
}


//Mitochondia; Intemediate Variable;  Optimize Things that use this; RT_over_F; exp_3_FRT_Dpsio;F_over_RT; 1/ADPm; ATPm
double calculateVATPase()
{
		double exp_3FRT_DmuH = exp( FRT_3 * DmuH );
		double AF1 =  kf1_Pi * ATPm / (ADPm);
		double denominator = - rhoF1 / ( exp_3_FRT_Dpsio + p1_exp_3_FRT_Dpsio * AF1 + ( p2 + p3 * AF1 ) * exp_3FRT_DmuH);

		VATPase = ( (VATPase_C1 + pc2 * exp_3FRT_DmuH) * AF1 - pa * exp_3FRT_DmuH ) * denominator;
		return VATPase;
}

//Mitochondia; Intemediate Variable;  Optimize Things that use this; RT_over_F; exp_3_FRT_Dpsio;F_over_RT; 1/ADPm; ATPm
double calculateVhu()
{
		double exp_3FRT_DmuH = exp( FRT_3 * DmuH );
		double AF1 =  kf1_Pi * ATPm / (ADPm);
		double denominator = - rhoF1 / ( exp_3_FRT_Dpsio + p1_exp_3_FRT_Dpsio * AF1 + ( p2 + p3 * AF1 ) * exp_3FRT_DmuH);

		Vhu     = ( pa_300 + pa_300 * AF1  - pa_pb_3 * exp_3FRT_DmuH ) * denominator;
		return Vhu;
}

double calculateVANT()
{
		//Parameters: gh, VmDT; hm
		double ATPi_ADP = ATPi / ADP;
		double ADPm_ATPm = (ADPm) / ATPm;

		VANT = ( VmDT_75 - VmDT_20 * ATPi_ADP * ADPm_ATPm * exp(-F_over_RT * (Dpsi)) ) /
				( ( 1.0 + tenDiv9 * ATPi_ADP * exp( -hm_F_over_RT * (Dpsi) ) ) * (1.0 + 18 * ADPm_ATPm ) );
		return VANT;
}


double calculateVhleak()
{
		//Parameters: gh, VmDT; hm
		Vhleak = gh * DmuH;
		return Vhleak;
}

// include VIMAC2, Vt2SO2m, VSOD, VCAT, VGPX, and VGR
	double calculateVIMAC2()
{
		VIMAC2 = -1.0 * (1.0e-3 + af / (1.0 + Kcc / (ROSi))) *
			(G_L + G_max / ( 1 + exp(Kappa * (Em - (-1.0 * (Dpsi)))))) * (-1.0 * (Dpsi));
		return VIMAC2;
}

	double calculateVt2SO2m()
{
		//unit of Vt2SO2m is mM/ms
		Vt2SO2m = -j_Vt2SO2m * ( (-1.0 * (Dpsi) - RT_over_F * log( (ROSm) / (ROSi))) / (Dpsi))
		* (-1 * (1.0e-3 + af / (1 + Kcc / (ROSi))) * (G_L + G_max / (1.0 + exp( Kappa * (Em - (-1.0 * (Dpsi)))))) * (-1.0 * (Dpsi)));
		return Vt2SO2m;
}

   
double calculateVSOD() {

    //new species ROSm (v12),ROSi (v13), H2O2 (v13),
    VSOD = 2 * k1_SOD * k5_SOD * (k1_SOD + k3_SOD * (1 + (H2O2 / Kki))) * etSOD * (ROSi) / (k5_SOD * (2 * k1_SOD + k3_SOD * (1 + (H2O2 / Kki))) + k1_SOD * k3_SOD * (1 + (H2O2 / Kki)) * (ROSi));

    return VSOD;
}


	double calculateVCAT()
{
		VCAT = 2 * k1_CAT * EtCAT * exp (-fr * (H2O2)) * (H2O2);
		return VCAT;
}

	double calculateVGPX()
{
		VGPX = EtGPX * (H2O2) * (GSH) / (Phi1 * (GSH) + Phi2 * (H2O2));//was the original one
		return VGPX;
}

	double calculateVGR()
{
		VGR = k1_GR * EtGR / (1 + Km_GSSG / (0.5 * (GT - (GSH))) + Km_NADPH / NADPH + Km_GSSG * Km_NADPH / (0.5 * (GT - (GSH)) / NADPH));//was the original one
		return VGR;
}

double calculateVuni(double Cai)
{
	//Parameters: Vmuni; ktrans; L; kact; na
		double Cai_ktrans_plus1 = 1.0 + (Cai) * inv_ktrans;
		double Cai_ktrans_plus1_p3 = Cai_ktrans_plus1 * Cai_ktrans_plus1 * Cai_ktrans_plus1;

		Vuni =  Vmuni_ktrans * (Cai) * FRT2_Dpsi * Cai_ktrans_plus1_p3  /
				( ( Cai_ktrans_plus1_p3 * Cai_ktrans_plus1 + L / pow( 1.0 + (Cai) * inv_kact, na) ) * (1.0 - exp(-FRT2_Dpsi)));
		return Vuni;
}

double calculateVnaca(double Cai)
{
	//Parameters: b; VmNC; Kna; n; Knca;
		Vnaca = VmNC * exp( b_05 * FRT2_Dpsi ) * (Cam) / ( (Cai) * pow( ( 1.0 + Kna / (Nai) ) , n) * ( 1.0 + Knca / (Cam) ) );
		return Vnaca;
}

//make void function bc it doesn't have return statement
void sharedVariableUpdate_Mitochondria() {
    //After calculateDmuH
    exp_FRT_6_g_DmuH = exp(FRT_6_g * DmuH);
    //Must be done after calculateNAD is Called
    KmIDNAD_NAD = KmIDNAD / NAD;
    //Anytime
    FRT2_Dpsi = FRT2 * ( (Dpsi) - 91.0 );
}
