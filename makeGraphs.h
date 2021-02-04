#ifndef MAKEGRAPHS_H
#define MAKEGRAPHS_H
void drawRectSB(double xmean, double ymean, double xstd, double ystd, int xsig, int xskip, int xsb, int ysig, int yskip, int ysb){
        // regions are:
        //  0 1 2
        //  3 4 5
        //  6 7 8
        //  where xmin, xmax, ymin,ymax all belong to region 5
        //  the 10th element is the interesection of the negation of all regions
        //  11th element is 1345 and 12th element is 0268
	TBox *box = new TBox();
	box->SetFillStyle(0);	
	box->SetLineWidth(2);
	box->SetLineColor(kRed);
	box->DrawBox(xmean-(xsig+xskip+xsb)*xstd,ymean+(ysig+yskip)*ystd,xmean-(xsig+xskip)*xstd,ymean+(ysig+yskip+ysb)*ystd); // 0
	box->DrawBox(xmean-xsig*xstd,ymean+(ysig+yskip)*ystd,xmean+xsig*xstd,ymean+(ysig+yskip+ysb)*ystd); // 1
	box->DrawBox(xmean+(xsig+xskip)*xstd,ymean+(ysig+yskip)*ystd,xmean+(xsig+xskip+xsb)*xstd,ymean+(ysig+yskip+ysb)*ystd); // 2
	box->DrawBox(xmean-(xsig+xskip+xsb)*xstd,ymean-ysig*ystd,xmean-(xsig+xskip)*xstd,ymean+ysig*ystd); // 3 
	box->DrawBox(xmean-xsig*xstd,ymean-ysig*ystd,xmean+xsig*xstd,ymean+ysig*ystd); // 4 
	box->DrawBox(xmean+(xsig+xskip)*xstd,ymean-ysig*ystd,xmean+(xsig+xskip+xsb)*xstd,ymean+ysig*ystd); // 5 
	box->DrawBox(xmean-(xsig+xskip+xsb)*xstd,ymean-(ysig+yskip+ysb)*ystd,xmean-(xsig+xskip)*xstd,ymean-(ysig+yskip)*ystd); // 6 
	box->DrawBox(xmean-xsig*xstd,ymean-(ysig+yskip+ysb)*ystd,xmean+xsig*xstd,ymean-(ysig+yskip)*ystd); // 7 
	box->DrawBox(xmean+(xsig+xskip)*xstd,ymean-(ysig+yskip+ysb)*ystd,xmean+(xsig+xskip+xsb)*xstd,ymean-(ysig+yskip)*ystd); // 8 
}


void drawLineRectSB(double xmin, double xmax, double xskip, double maxValue){
	double xlength = xmax-xmin;
	TLine *line = new TLine( xmin, 0, xmin, maxValue);
	line->Draw();
	line->SetLineColor(kMagenta);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->DrawLine( xmax, 0, xmax, maxValue );
	line->DrawLine( xmax+xskip, 0, xmax+xskip, maxValue );
	line->DrawLine( xmax+xskip+xlength/2, 0, xmax+xskip+xlength/2, maxValue );
	line->DrawLine( xmin-xskip, 0, xmin-xskip, maxValue );
	line->DrawLine( xmin-xskip-xlength/2, 0, xmin-xskip-xlength/2, maxValue );
}

#endif 
