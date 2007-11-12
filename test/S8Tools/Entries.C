
#include<iomanip>

void Entries() {

	int entries = 0;
	int total = 0;
	int entries2 =0;
	//cout << std::fixed;
	
	TChain f1("summary");
	f1.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_BBbar/lpc_BBbar_30_50_*.root");
	entries = f1.GetEntries();	
	cout << "semileptonic \\bb 30-50  & " << entries << "\\\\" << endl;
	f1.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_BBbar/lpc_BBbar_50_80_*.root");
	entries2 = f1.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\bb 50-80  & " << entries2 << "\\\\" << endl;
	f1.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_BBbar/lpc_BBbar_80_120_*.root");
	entries2 = f1.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\bb 80-120  & " << entries2 << "\\\\" << endl;
	f1.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_BBbar/lpc_BBbar_120_170_*.root");
	entries2 = f1.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\bb 120-170  & " << entries2 << "\\\\" << endl;
	f1.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_BBbar/lpc_BBbar_170_230_*.root");
	entries2 = f1.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\bb 170-230  & " << entries2 << "\\\\" << endl;
	
	total = f1.GetEntries();

	TChain f2("summary");
	f2.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_CCbar/lpc_CCbar_30_50_*.root");
	entries = f2.GetEntries();
	cout << "semileptonic \\cc 30-50  & " << entries << "\\\\" << endl;
	f2.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_CCbar/lpc_CCbar_50_80_*.root");
	entries2 = f2.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\cc 50-80  & " << entries2 << "\\\\" << endl;
	f2.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_CCbar/lpc_CCbar_80_120_*.root");
	entries2 = f2.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\cc 80-120  & " << entries2 << "\\\\" << endl;
	f2.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_CCbar/lpc_CCbar_120_170_*.root");
	entries2 = f2.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\cc 120-170  & " << entries2 << "\\\\" << endl;
	f2.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/lpc_CCbar/lpc_CCbar_170_230_*.root");
	entries2 = f2.GetEntries() - entries;
	entries += entries2;
	cout << "semileptonic \\cc 170-230  & " << entries2 << "\\\\" << endl;
	
	total += f2.GetEntries();

	TChain f3("summary");
	f3.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/QCD_20_30/*.root");
	entries = f3.GetEntries();
	cout << "QCD $\hat{p_T}$ 20-30  & " << entries << "\\\\" << endl;
	f3.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/QCD_30_50/*.root");
	entries2 = f3.GetEntries() - entries;
	entries += entries2;
	cout << "QCD $\hat{p_T}$ 30-50  & " << entries2 << "\\\\" << endl;
	f3.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/QCD_50_80/*.root");
	entries2 = f3.GetEntries() - entries;
	entries += entries2;
	cout << "QCD $\hat{p_T}$ 50-80  & " << entries2 << "\\\\" << endl;
	f3.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/QCD_80_120/*.root");
	entries2 = f3.GetEntries() - entries;
	entries += entries2;
	cout << "QCD $\hat{p_T}$ 80-120  & " << entries2 << "\\\\" << endl;
	f3.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/QCD_120_170/*.root");
	entries2 = f3.GetEntries() - entries;
	entries += entries2;
	cout << "QCD $\hat{p_T}$ 120-170  & " << entries2 << "\\\\" << endl;
	f3.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/QCD_170_230/*.root");
	entries2 = f3.GetEntries() - entries;
	entries += entries2;
	cout << "QCD $\hat{p_T}$ 170-230  & " << entries2 << "\\\\" << endl;
	f3.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/QCD_230_300/*.root");
	entries2 = f3.GetEntries() - entries;
	entries += entries2;
	cout << "QCD $\hat{p_T}$ 230-300  & " << entries2 << "\\\\" << endl;
	total += f3.GetEntries();
	
	cout << "Soup (Sum of above)  & " << total << "\\\\" << endl;
	
	total = 0;
	
	TChain f4("summary");
	f4.Add("/uscms/home/jandrea/OctDataSys8/InclBBar_30_50/*.root");
	entries = f4.GetEntries();
	cout << "Incl \\bb 30-50  & " << entries << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/InclBBar_50_80/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "Incl \\bb 50-80  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/InclBBar_80_120/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "Incl \\bb 80-120  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/InclBBar_120_170/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "Incl \\bb 120-170  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/InclBBar_170up/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "Incl \\bb 170-up  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/BBbar_50_80/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "inclusive \\bb 50-80  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/BBbar_80_120/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "inclusive \\bb 90-120  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/CCbar_30_50/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "inclusive \\cc 30-50  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/CCbar_50_80/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "inclusive \\cc 50-80  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/CCbar_80_120/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "inclusive \\cc 80-120  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/CCbar_120_170/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "inclusive \\cc 120-170  & " << entries2 << "\\\\" << endl;
	f4.Add("/uscms/home/jandrea/OctDataSys8/CCbar_170_230/*.root");
	entries2 = f4.GetEntries() - entries;
	entries += entries2;
	cout << "inclusive \\cc 170-230  & " << entries2 << "\\\\" << endl;

	total = f4.GetEntries();

	TChain f5("summary");
	f5.Add("~/lpcbtag/PerformanceMeasurements/2007_Oct_9_13X/TTbar_inclusive_onemuon_TopRex/*.root");
	entries = f5.GetEntries();
	cout << "inclusive \ttbar & " << entries << "\\\\" << endl;
}
