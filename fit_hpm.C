#include "TFile.h"
#include "TMath.h"
#include "TH3F.h"
#include "TF3.h"
#include "TH2D.h"
#include "TF2.h"

double gaus3d(double *X, double *a) {
    double N = a[0];
    double m_x = a[1];
    double m_y = a[2];
    double m_z = a[3];

    double K_xx = a[4];
    double K_yy = a[5];
    double K_zz = a[6];

    double K_xy = a[7];
    double K_yz = a[8];
    double K_zx = a[9];

    double x = X[0] - m_x;
    double y = X[1] - m_y;
    double z = X[2] - m_z;

    return N * TMath::Exp(
        -(x*x*K_xx + y*y*K_yy + z*z*K_zz + 2*x*y*K_xy + 2*y*z*K_yz + 2*z*x*K_zx)
    );

}

double gaus2d(double *X, double *a) {
    double N = a[0];
    double m_x = a[1];
    double m_y = a[2];

    double K_xx = a[3];
    double K_yy = a[4];

    double K_xy = a[5];

    double x = X[0] - m_x;
    double y = X[1] - m_y;

    return N * TMath::Exp(
        -(x*x*K_xx + y*y*K_yy + 2*x*y*K_xy)
    );

}

void fit_hpm() {
    TFile *file = new TFile("pm_0.root","read");
    TH3F *h = (TH3F*)file->Get("hpm");
/*
    TF3 *f3 = new TF3("f3", gaus3d, 0.0, 1.0, 0.0, 2*TMath::Pi(), 0.0, 2.0, 10);

    h->Fit(f3);
    */
    TF2 *f2 = new TF2("f2", gaus2d, 0.4, 1.4, 4.0, 6.0, 6);

    h->GetXaxis()->SetRange(86,86);
    TH2D *p_yz = (TH2D*)h->Project3D("yz");
    p_yz->Draw("CONT3");

    f2->SetParameter(1, 0.8);
    f2->SetParameter(2, 5.0);

    p_yz->Fit(f2,"R");
}