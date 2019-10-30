#ifndef PLOT_H
#define PLOT_H

#include <assert.h>

class loadAmplitudes{
	private:
		int _numAmps=8;
		int _numBins=65;
		ifstream file;
		// seems like we cant reserve in a class without being in a function to be called later
	public:
		std::vector<double> masses;
		std::vector<double> S0ms;
		std::vector<double> P0ms;
		std::vector<double> P1ms;
		std::vector<double> D0ms;
		std::vector<double> D1ms;
		std::vector<double> P1ps;
		std::vector<double> D1ps;
		std::vector<double> totals;
		map < string, std::vector<double> > m_nameToVect;

        	loadAmplitudes ( string fileLoc, int numAmps, int numBins){
			cout << "Loaded " << fileLoc << endl;
			file.open(fileLoc);	
			assert ( _numAmps == numAmps );	
			assert ( _numBins == numBins );	
		}
		double mass;
		double S0m[2], P0m[2], P1m[2], D0m[2], D1m[2], P1p[2], D1p[2], total[2];

		void load(){
			S0ms.reserve(_numBins*2);
			P0ms.reserve(_numBins*2);
			P1ms.reserve(_numBins*2);
			D0ms.reserve(_numBins*2);
			D1ms.reserve(_numBins*2);
			P1ps.reserve(_numBins*2);
			D1ps.reserve(_numBins*2);
			totals.reserve(_numBins*2);
			while ( file >> mass >> S0m[0] >> S0m[1] >> P0m[0] >> P0m[1] >> P1m[0] >> P1m[1] >> D0m[0] >> D0m[1]>> D1m[0] >> D1m[1]>> P1p[0] >> P1p[1]>> D1p[0] >> D1p[1] 
					>> total[0] >> total[2]) {
				masses.push_back(mass);
				S0ms.push_back(S0m[0]);
				P0ms.push_back(P0m[0]);
				P1ms.push_back(P1m[0]);
				D0ms.push_back(D0m[0]);
				D1ms.push_back(D1m[0]);
				P1ps.push_back(P1p[0]);
				D1ps.push_back(D1p[0]);
				totals.push_back(total[0]); 
				S0ms.push_back(S0m[1]);
				P0ms.push_back(P0m[1]);
				P1ms.push_back(P1m[1]);
				D0ms.push_back(D0m[1]);
				D1ms.push_back(D1m[1]);
				P1ps.push_back(P1p[1]);
				D1ps.push_back(D1p[1]);
				totals.push_back(total[1]); 
			}
			m_nameToVect["S0m"] = S0ms; 
			m_nameToVect["P0m"] = P0ms; 
			m_nameToVect["P1m"] = P1ms; 
			m_nameToVect["D0m"] = D0ms; 
			m_nameToVect["D1m"] = D1ms; 
			m_nameToVect["P1p"] = P1ps; 
			m_nameToVect["D1p"] = D1ps; 
			m_nameToVect["total"] = totals; 
		}

		std::vector<double> getAmp(string nameAmp){
			return m_nameToVect[nameAmp];
		}

		
};

#endif
