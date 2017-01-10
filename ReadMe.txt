Possible Errors:

(1) error LNK1104: cannot open file '..\mergeBathy_CPP\x64\Release\mergeBathy.exe' 

	Build > Clean mergeBathy



(2) error LNK1123: failure during conversion to COFF: file invalid or corrupt

	The solution was found via googling at:
	http://stackoverflow.com/questions/10888391/error-link-fatal-error-lnk1123-failure-during-conversion-to-coff-file-inval

	In summary, you can either:
	(A) Reintstall VS 2010 SP1 and/or VS 2010 Compiler Pack

	(B) or the simplest, make sure your using cvtres.exe in the first location below and not the second.  This requires either adjusting your path or renaming the file in the second location. 
	"You may have two versions of the cvtres.exe utility in your path. 
	One at C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\BIN\cvtres.exe 
	and one at C:\Windows\Microsoft.NET\Framework\v4.0.30319\cvtres.exe"



