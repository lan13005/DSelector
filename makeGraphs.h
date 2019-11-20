#ifndef MAKEGRAPHS_H
#define MAKEGRAPHS_H

void drawRectSB(double xmin, double xmax, double ymin, double ymax, double xskip, double yskip){
        // regions are:
        //  0 1 2
        //  3 4 5
        //  6 7 8
        //  where xmin, xmax, ymin,ymax all belong to region 5
        //  the 10th element is the interesection of the negation of all regions
        //  11th element is 1345 and 12th element is 0268
        double xlength = xmax-xmin;
        double ylength = ymax-ymin;
	TBox *box = new TBox(xmin,ymin,xmax,ymax);
	box->SetFillStyle(0);	
	box->SetLineColor(kRed);
	box->Draw(); // 4 
	box->DrawBox(xmax+xskip,ymin,xmax+xskip+xlength/2,ymax); // 5
	box->DrawBox(xmin-xskip-xlength/2,ymin,xmin-xskip,ymax); // 3
	box->DrawBox(xmin,ymax+yskip,xmax,ymax+yskip+ylength); // 1
	box->DrawBox(xmax+xskip,ymax+yskip,xmax+xskip+xlength/2,ymax+yskip+ylength); //2
	box->DrawBox(xmin-xskip-xlength/2,ymax+yskip,xmin-xskip,ymax+yskip+ylength); //0
	box->DrawBox(xmin,ymin-yskip-ylength,xmax,ymin-yskip); // 7
	box->DrawBox(xmin-xskip-xlength/2,ymin-yskip-ylength,xmin-xskip,ymin-yskip); // 6
	box->DrawBox(xmax+xskip,ymin-yskip-ylength,xmax+xskip+xlength/2,ymin-yskip); // 8
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
