package alignment;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.InputMismatchException;
import java.util.Map;
import java.util.Scanner;

public class SequenceAlignment {
	
	private int[][] calcTable;			// nasza glowna tabela
	private int[][][] indicatorsTable;	//przechowuje wartosci indeksow skad pochodza wartosci w calcTable

	private Map<Character, Integer> gapPenalty = new HashMap<Character, Integer>(){}; // mapa przechowujaca wartosci kar za indel pobrane z macierzy
	static final String GAP_CHAR = "_"; //znak indela
	private static Map<Integer, String> MemRNASeq1 = new HashMap<Integer, String>(){}; // przechowuje trojki RNA ktore zostaja zamienione na aminokwasy aby ponownie przeksztalcic je w RNA
	private static Map<Integer, String> MemRNASeq2 = new HashMap<Integer, String>(){}; // przechowuje trojki RNA ktore zostaja zamienione na aminokwasy aby ponownie przeksztalcic je w RNA
	
	public SequenceAlignment(Map<Character, Integer> gapPenalty) {
		this.gapPenalty = gapPenalty;
	}


	public void calcGlobalAlignment(String sequence1Original, String sequence2Original, boolean printResults, int [][]matrix, int matrixType, int isRNA) {
		String seq1 = normalizeSeq(sequence1Original);
		String seq2 = normalizeSeq(sequence2Original);

		//Deklaracja rozmiaru tabel
		calcTable = new int[seq2.length()][seq1.length()];
		indicatorsTable = new int[seq2.length()][seq1.length()][2];
		
		//wypelnienie miejsca 0,0 
		calcTable[0][0] = 0; 
		
		//wypelnienie 0 kolumny
		for (int i = 1; i < seq2.length(); i++) {
			try {
				calcTable[i][0] = i * this.gapPenalty.get(seq2.charAt(i));
			}
			catch (NullPointerException e) {
				System.out.println("\nBlad podczas wykorzystania wartosci z macierzy\n");
				System.exit(404);
			}
			indicatorsTable[i][0][0] = i - 1;
			indicatorsTable[i][0][1] = 0;
		}
		//wypelnienie 0 wiersza
		for (int j = 1; j < seq1.length(); j++) {
			try {
				calcTable[0][j] = j * this.gapPenalty.get(seq1.charAt(j));
			}
			catch (NullPointerException e) {
				System.out.println("\nBlad podczas wykorzystania wartosci z macierzy\n");
				System.exit(404);
			}
			indicatorsTable[0][j][0] = 0;
			indicatorsTable[0][j][1] = j - 1;
		}
		// brak wskaznika w rogu tabeli
		indicatorsTable[0][0][0] = -1;
		indicatorsTable[0][0][1] = -1;


		//wypelnienie wnetrza tabeli
		for (int j = 1; j < seq1.length(); j++) {
			for (int i = 1; i < seq2.length(); i++) {
				//System.out.println(i+" "+j);	
				int normalAlign = matrix[i-1][j-1] + calcTable[i - 1][j - 1];	//case1: seq1 i seq2 normalnie dopasowane
				int seq1WithGap = -9999;
				int seq2WithGap = -9999;
				try {
					seq1WithGap = this.gapPenalty.get(seq1.charAt(j)) + calcTable[i][j - 1];		//case2: seq1 z indelem
					seq2WithGap = this.gapPenalty.get(seq2.charAt(i)) + calcTable[i - 1][j];		//case3: seq2 z indelem
		
				} catch (NullPointerException e) {
					System.out.println("\nBlad podczas wykorzystania wartosci z macierzy\n");
					System.exit(404);
				}
				
				if(matrixType == 1) { // macierz odleglosci
					//obliczanie min
					if (normalAlign <= seq1WithGap && normalAlign <= seq2WithGap) {			//case1 jest min
						calcTable[i][j] = normalAlign;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j - 1;
					}
					else if (seq2WithGap <= normalAlign && seq2WithGap <= seq1WithGap) {	//case3 jest min
						calcTable[i][j] = seq2WithGap;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j;
					}
					else {																	//case2 jest min
						calcTable[i][j] = seq1WithGap;
						indicatorsTable[i][j][0] = i;
						indicatorsTable[i][j][1] = j - 1;
					}
				}
				
				else if(matrixType == 2) { // macierz podobienstwa
					//obliczanie max
					if (normalAlign >= seq1WithGap && normalAlign >= seq2WithGap) {//case1 jest max
						calcTable[i][j] = normalAlign;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j - 1;
					}
					else if (seq2WithGap >= normalAlign && seq2WithGap >= seq1WithGap) {	//case3 jest max
						calcTable[i][j] = seq2WithGap;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j;
					}
					else {																	//case2 jest max
						calcTable[i][j] = seq1WithGap;
						indicatorsTable[i][j][0] = i;
						indicatorsTable[i][j][1] = j - 1;
					}
				}
			}
		}

		if(printResults){
			System.out.println("\nZestawienie \""+sequence1Original+"\" z \""+sequence2Original+"\"");
			System.out.println("\nTablica\n");
			print2DimTable(calcTable);
			System.out.println("\nTablica wskzanikow\n");	
			printIndicatorsTable(indicatorsTable);

			int minimumPenalty = calcTable[sequence2Original.length()][sequence1Original.length()];
	 		System.out.println("\n" + minimumPenalty + "     jest maksymalna wartoscia dopasowania \""+sequence1Original+"\" z \""+sequence2Original+"\"");
			findGlobalAlignment(seq1, seq2,isRNA);
		}
	}
	
	public void calcLocalAlignment(String sequence1Original, String sequence2Original, boolean printResults, int [][]matrix, int matrixType,int isRNA) {
		String seq1 = normalizeSeq(sequence1Original);
		String seq2 = normalizeSeq(sequence2Original);

		//Deklaracja rozmiaru tabel
		calcTable = new int[seq2.length()][seq1.length()];
		indicatorsTable = new int[seq2.length()][seq1.length()][2];

		//wypelnienie 0 kolumny
		for (int i = 0; i < seq2.length(); i++) {	
			calcTable[i][0] = 0;
		}
		//wypelnienie 0 wiersza
		for (int j = 0; j < seq1.length(); j++) {	
			calcTable[0][j] = 0;
		}
		//Na poczatu ustawiam cala tablice wskaznikow na -1 czyli brak wskaznika
		for (int i = 0; i < seq2.length(); i++) {	
			for (int j = 0; j < seq1.length(); j++) {	
				indicatorsTable[i][j][0] = -1;
				indicatorsTable[i][j][1] = -1;
			}
		}


		//wypelnienie wnetrza tabelie
		for (int j = 1; j < seq1.length(); j++) {
			for (int i = 1; i < seq2.length(); i++) {
				//System.out.println(i+" "+j);	
				int normalAlign = matrix[i-1][j-1] + calcTable[i - 1][j - 1];	//case1: seq1 i seq2 normalnie dopasowane
				int seq1WithGap = -9999;
				int seq2WithGap = -9999;
				try {
					seq1WithGap = this.gapPenalty.get(seq1.charAt(j)) + calcTable[i][j - 1];		//case2: seq1 z indelem
					seq2WithGap = this.gapPenalty.get(seq2.charAt(i)) + calcTable[i - 1][j];		//case3: seq2 z indelem
		
				} catch (NullPointerException e) {
					System.out.println("\nBlad podczas wykorzystania wartosci z macierzy\n");
					System.exit(404);
				}
				
				int zero = 0;													//case4: 0
				
				if(matrixType == 1) { // macierz odleglosci
					//obliczanie minimum
					if (normalAlign <= seq1WithGap && normalAlign <= seq2WithGap && normalAlign <= 0) {//case1 jest min
						calcTable[i][j] = normalAlign;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j - 1;
					}
					else if (seq2WithGap <= normalAlign && seq2WithGap <= seq1WithGap && seq2WithGap <= 0) {	//case3 jest min
						calcTable[i][j] = seq2WithGap;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j;
					}
					else if (seq1WithGap <= normalAlign && seq1WithGap <= seq2WithGap && seq1WithGap <= 0){	//case2 jest min
						calcTable[i][j] = seq1WithGap;
						indicatorsTable[i][j][0] = i;
						indicatorsTable[i][j][1] = j - 1;
					}
					else { 																					// case4 jest min
						calcTable[i][j]= 0;  // jesli min jest 0 to nie ustawiam wskaznika, czyli zostaje -1,-1
					}
				}
				
				else if(matrixType == 2) { // macierz podobienstwa
					//obliczanie maximum
					if (normalAlign >= seq1WithGap && normalAlign >= seq2WithGap && normalAlign >= 0) {//case1 jest max
						calcTable[i][j] = normalAlign;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j - 1;
					}
					else if (seq2WithGap >= normalAlign && seq2WithGap >= seq1WithGap && seq2WithGap >= 0) {	//case3 jest max
						calcTable[i][j] = seq2WithGap;
						indicatorsTable[i][j][0] = i - 1;
						indicatorsTable[i][j][1] = j;
					}
					else if (seq1WithGap >= normalAlign && seq1WithGap >= seq2WithGap && seq1WithGap >= 0){	//case2 jest max
						calcTable[i][j] = seq1WithGap;
						indicatorsTable[i][j][0] = i;
						indicatorsTable[i][j][1] = j - 1;
					}
					else { 																					// case4 jest max
						calcTable[i][j]= 0;   // jesli max jest 0 to nie ustawiam wskaznika, czyli zostaje -1,-1
					}
				}
			}
		}

		if(printResults){
			System.out.println("\nZestawienie \""+sequence1Original+"\" z \""+sequence2Original+"\"");
			System.out.println("\nTablica\n");
			print2DimTable(calcTable);
			System.out.println("\nTablica wskzanikow\n");	
			printIndicatorsTable(indicatorsTable);
			
			int [] maxValueIndexArray = new int [2];
			int minimumPenalty = findMaximumIn2DArrayWithIndex(calcTable, maxValueIndexArray);  // szukam wartosci najwieszkej w tabeli i przechowuje jej indeks
	 		System.out.println("\n" + minimumPenalty + "   jest maksymalna wartoscia dopasowania \""+sequence1Original+"\" z \""+sequence2Original+"\"");
	 		//System.out.println("Max value index" + maxValueIndexArray[0] + "   " + maxValueIndexArray[1]);	
	 		findLocalAlignment(seq1, seq2,maxValueIndexArray,isRNA);
		}
	}
	
	//znajduje najwiekszy element w tablicy 2D i zapisuje jego indeks
	private int findMaximumIn2DArrayWithIndex(int[ ][ ] a, int []answerArray)
	{
	    int maxVal = -99999;
	    for(int row = 0; row < a.length; row++)
	    {
	        for(int col = 0; col < a[row].length; col++)
	        {
	            if(a[row][col] > maxVal)
	            {
	                maxVal = a[row][col];
	                answerArray[0] = row;
	                answerArray[1] = col;
	            }
	        }
	    }
	    return maxVal;
	}

	private void print2DimTable(int[][] table) {
		for (int[] row : table) {
			for (int value : row) {
				System.out.print(value + "\t");
			}
			System.out.println();
		}
	}

	private void printIndicatorsTable(int[][][] table3D) {
		for (int[][] row : table3D) {
			for (int[] xyPair : row) {
				System.out.print(Arrays.toString(xyPair) + "  \t");
			}
			System.out.println();
		}
	}

 // znajduje i wypisuje optymalne globalne dopasowanie 
	private void findGlobalAlignment(String seq1, String seq2, int isRNA) {
		String seq1Aligned = ""; 	
		String seq2Aligned = "";
		
		//iteratory po wierszach i kolumnach tabeli
		int i = seq2.length() - 1; //-1 zeby pominac bialy znak z poczatku wyrazu
		int j = seq1.length() - 1;
		
		//zmienne pomocnicze do przechodzenia po wyrazach
		int seq1L = seq1.length() - 1; 
		int seq2L = seq2.length() - 1;
		
		String seq1RNAFromAmino = ""; // wykorzystywane w przypadku zamiany aminokwasow na RNA
		String seq2RNAFromAmino = "";
		
		while (i > 0 || j > 0) {
			//System.out.println(predecessorIndexes[i][j][0]+"   "+predecessorIndexes[i][j][1]);
			
			//indel w seq1
			if(indicatorsTable[i][j][0] == i-1 && indicatorsTable[i][j][1] == j) {
				
				seq2Aligned = seq2.charAt(seq2L) + seq2Aligned;
				seq1Aligned = GAP_CHAR + seq1Aligned;
				seq2RNAFromAmino = MemRNASeq2.get(seq2L-1) + seq2RNAFromAmino;
				seq1RNAFromAmino = "___" + seq1RNAFromAmino;
				i--;
				seq2L--;
				
			}
			//indel w seq2
			else if(indicatorsTable[i][j][0] == i && indicatorsTable[i][j][1] == j-1) {
				
				seq1Aligned = seq1.charAt(seq1L) + seq1Aligned;
				seq2Aligned = GAP_CHAR + seq2Aligned;
				seq1RNAFromAmino = MemRNASeq1.get(seq1L-1) + seq1RNAFromAmino;
				seq2RNAFromAmino = "___" + seq2RNAFromAmino;
				j--;
				seq1L--;
			}
			//dopasowanie po przek¹tnej
			else if(indicatorsTable[i][j][0] == i-1 && indicatorsTable[i][j][1] == j-1) {
				
				seq1Aligned = seq1.charAt(seq1L) + seq1Aligned;
				seq2Aligned = seq2.charAt(seq2L) + seq2Aligned;
				seq1RNAFromAmino = MemRNASeq1.get(seq1L-1) + seq1RNAFromAmino;
				seq2RNAFromAmino = MemRNASeq2.get(seq2L-1) + seq2RNAFromAmino;
				i--;
				j--;
				seq1L--;
				seq2L--;
			}		

			//System.out.println(i+"   "+j+"\n" + seq1Aligned + "\n" + seq2Aligned + "\n");

		}
		if(isRNA == 1)
			System.out.println("\nOptymalne dopasowanie globalne:\n\n" + seq1Aligned + "\n" + seq2Aligned + "\n");
		else if (isRNA ==2)
		{
			System.out.println("\nOptymalne dopasowanie globalne w alfabecie aminokwasowym:\n\n" + seq1Aligned + "\n" + seq2Aligned + "\n");
			
			System.out.println("\nOptymalne dopasowanie globalne w alfabecie RNA:\n\n" + seq1RNAFromAmino + "\n" + seq2RNAFromAmino + "\n");
		}
		
	}
	 // znajduje i wypisuje optymalne lokalne dopasowanie
	private void findLocalAlignment(String seq1, String seq2, int [] maxIndex, int isRNA) {
		String seq1Aligned = ""; 	//Holds the actual sequence with gaps added
		String seq2Aligned = "";

		//iteratory po wierszach i kolumnach tabeli
		int i = maxIndex[0]; 	//-1 zeby pominac bialy znak z poczatku wyrazu
		int j = maxIndex[1];
		
		//zmienne pomocnicze do przechodzenia po wyrazach
		int seq1L = maxIndex[1]; 
		int seq2L = maxIndex[0];
		
		String seq1RNAFromAmino = ""; // wykorzystywane w przypadku zamiany aminokwasow na RNA
		String seq2RNAFromAmino = "";
		
		while (i > 0 || j > 0) {
			//System.out.println(predecessorIndexes[i][j][0]+"   "+predecessorIndexes[i][j][1]);
			
			if(indicatorsTable[i][j][0] != -1 && indicatorsTable[i][j][1] != -1) { // dopoki nie osiagnieto pola bez wskaznika
				
				//indel w seq1
				if(indicatorsTable[i][j][0] == i-1 && indicatorsTable[i][j][1] == j) {
					
					seq2Aligned = seq2.charAt(seq2L) + seq2Aligned;
					seq1Aligned = GAP_CHAR + seq1Aligned;
					seq2RNAFromAmino = MemRNASeq2.get(seq2L-1) + seq2RNAFromAmino;
					seq1RNAFromAmino = "___" + seq1RNAFromAmino;
					i--;
					seq2L--;
				}
				//indel w seq2
				else if(indicatorsTable[i][j][0] == i && indicatorsTable[i][j][1] == j-1) {
					
					seq1Aligned = seq1.charAt(seq1L) + seq1Aligned;
					seq2Aligned = GAP_CHAR + seq2Aligned;
					seq1RNAFromAmino = MemRNASeq1.get(seq1L-1) + seq1RNAFromAmino;
					seq2RNAFromAmino = "___" + seq2RNAFromAmino;
					j--;
					seq1L--;
				}
				//dopasowanie po przek¹tnej
				else if(indicatorsTable[i][j][0] == i-1 && indicatorsTable[i][j][1] == j-1) {
					
					seq1Aligned = seq1.charAt(seq1L) + seq1Aligned;
					seq2Aligned = seq2.charAt(seq2L) + seq2Aligned;
					seq1RNAFromAmino = MemRNASeq1.get(seq1L-1) + seq1RNAFromAmino;
					seq2RNAFromAmino = MemRNASeq2.get(seq2L-1) + seq2RNAFromAmino;
					i--;
					j--;
					seq1L--;
					seq2L--;
				}		
			}
			
			else { // osiagnieto pole bez wskaznika
				break;
			}
			//System.out.println(i+"   "+j+"\n" + seq1Aligned + "\n" + seq2Aligned + "\n");

		}
		
		if(isRNA == 1)
			System.out.println("\nOptymalne dopasowanie lokalne:\n\n" + seq1Aligned + "\n" + seq2Aligned + "\n");
		else if (isRNA ==2)
		{
			System.out.println("\nOptymalne dopasowanie lokalne w alfabecie aminokwasowym:\n\n" + seq1Aligned + "\n" + seq2Aligned + "\n");
			
			System.out.println("\nOptymalne dopasowanie lokalne w alfabecie RNA:\n\n" + seq1RNAFromAmino + "\n" + seq2RNAFromAmino + "\n");
		}
	}


	private String normalizeSeq(String sequence){
		return " " + sequence		//pusty znak na poczatku aby dodac 1 wiersz i kolumne do tablicy. 
				.toUpperCase();
	}
	
	public static String translateRNA(String seq, int whichSeq) {
		
	    String seq1 = seq.toUpperCase();
	    String result = "";
	    
		Map<String, Character> dictionary = new HashMap<String, Character>(){{
	        put("UUU", 'F');put("UCU", 'S');put("UAU", 'Y');put("UGU", 'C');
	        put("UUC", 'F');put("UCC", 'S');put("UAC", 'Y');put("UGC", 'C');
	        put("UUA", 'L');put("UCA", 'S');put("UAA", '-');put("UGA", '-');
	        put("UUG", 'L');put("UCG", 'S');put("UAG", '-');put("UGG", 'W');
	        put("CUU", 'L');put("CCU", 'P');put("CAU", 'H');put("CGU", 'R');
	        put("CUC", 'L');put("CCC", 'P');put("CAC", 'H');put("CGC", 'R');
	        put("CUA", 'L');put("CCA", 'P');put("CAA", 'Q');put("CGA", 'R');
	        put("CUG", 'L');put("CCG", 'P');put("CAG", 'Q');put("CGG", 'R');
	        put("AUU", 'I');put("ACU", 'T');put("AAU", 'N');put("AGU", 'S');
	        put("AUC", 'I');put("ACC", 'T');put("AAC", 'N');put("AGC", 'S');
	        put("AUA", 'I');put("ACA", 'T');put("AAA", 'K');put("AGA", 'R');
	        put("AUG", 'M');put("ACG", 'T');put("AAG", 'K');put("AGG", 'R');
	        put("GUU", 'V');put("GCU", 'A');put("GAU", 'D');put("GGU", 'G');
	        put("GUC", 'V');put("GCC", 'A');put("GAC", 'D');put("GGC", 'G');
	        put("GUA", 'V');put("GCA", 'A');put("GAA", 'E');put("GGA", 'G');
	        put("GUG", 'V');put("GCG", 'A');put("GAG", 'E');put("GGG", 'G');
	    }};

	    //zamiana trojek nukleotydowych na aminokwas na podstawie mapy
	    for (int i=0; i < Math.floor(seq1.length() / 3) * 3; i+=3)
	    {
    		String pomString = String.valueOf(seq1.charAt(i)) + String.valueOf(seq1.charAt(i+1)) + String.valueOf(seq1.charAt(i+2));
    				
	    	if(dictionary.containsKey(pomString)) {
	    		
			    char a = dictionary.get(pomString);
			    result += String.valueOf(a);
			    
			    if(whichSeq == 1)
			    	MemRNASeq1.put(i/3,pomString);
			    else if(whichSeq == 2)
			    	MemRNASeq2.put(i/3,pomString);
				//System.out.println(result);			    
	    	}
	    	else {
	    		System.out.println("\nBlad podczas translacji\n");
				System.exit(404);
	    	}
	    }
	    //System.out.println(MemRNASeq1+"\n"+MemRNASeq2);	

	    
	    return result;
	}
	
	
	
	public static void main(String[] args) {
		
		//czytanie 1 sekwencji z pliku
		Scanner sc1 = null;
		System.out.println("Czytam sekwencje 1 z pliku");
		try {
			sc1 = new Scanner(new File(".\\src\\samples\\seqLoc11.txt")); // nazwa pliku z 1 sekwencja
		} catch (FileNotFoundException e) {

			e.printStackTrace();
			System.out.println("\nBlad podczas odczytu\n");
			System.exit(404);
		} 
			    
		//czytanie 2 sekwencji z pliku
		Scanner sc2 = null;
		System.out.println("Czytam sekwencje 2 z pliku");
		try {
			sc2 = new Scanner(new File(".\\src\\samples\\seqLoc22.txt")); // nazwa pliku z 2 sekwencja
		} catch (FileNotFoundException e) {

			e.printStackTrace();
			System.out.println("\nBlad podczas odczytu\n");
			System.exit(404);
		}   
		
		//czytanie macierzy z pliku
		Scanner sc3 = null;
		System.out.println("Czytam macierz z pliku");
		try {
			sc3 = new Scanner(new File(".\\src\\samples\\matrixLocSim1122.txt")); // nazwa pliku z macierza
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.out.println("\nBlad podczas odczytu\n");
			System.exit(404);
		}
		
		int rows = sc3.nextInt();			//czytanie rozmiarów macierzy
		int columns = sc3.nextInt();
		
		int [][] matrixRead = new int[rows][columns];

		char[] seq1Pom = new char [columns]; 	// sekwencje pobrane z pliku z macierz¹
		char [] seq2Pom = new char [rows];
		
		sc3.nextLine();
		//sc3.next();
		for(int j = 0; j < columns; ++j)
	    {
			seq1Pom[j] = sc3.next().toUpperCase().charAt(0);
	    }
		sc3.nextLine();
		
		//czytanie macierzy
		for(int i = 0; i < rows; ++i)
		{
			seq2Pom[i] = sc3.next().toUpperCase().charAt(0);
		    for(int j = 0; j < columns; ++j)
		    {
		        if(sc3.hasNextInt())
		        {
		            matrixRead[i][j] = sc3.nextInt();
		        }
		    }
		}
		String seq1 = sc1.nextLine().toUpperCase(); // przypisanie sekwencji pobranych z plików do zmiennych
		String seq2 = sc2.nextLine().toUpperCase();
		
		System.out.println("\nSekwencja pierwsza:\t"+seq1);
		System.out.println("Sekwencja druga:\t"+seq2);
		System.out.println("\nMacierz: \n");
		System.out.print("\t");
		
		// wypisanie macierzy
		for(int i = 0; i < columns; ++i)
			System.out.print(seq1Pom[i]+"\t");
		System.out.println("\n");
		
		for(int i = 0; i < rows; ++i)
		{
			System.out.print(seq2Pom[i]+"\t");
			
		    for(int j = 0; j < columns; ++j)
		    {
		        
		    		System.out.print(matrixRead[i][j]+"\t");
		        
		    }
		    System.out.println();
		}
		
		
		//wprowadzanie typu dopasowaniua z linii polecen
		Scanner scan = null;
				boolean isSuccess4 = false;
				int isRNA = 0;
				System.out.println("\nDopasowac sekwencje jako RNA ?: 1 - nie, 2 - tak\n");
				scan = null;
				do 
				{
					try {
						scan = new Scanner(System.in);
						isRNA = scan.nextInt();
					} catch (InputMismatchException e) {
						isSuccess4 = false;
					}
					
					if(isRNA == 1 || isRNA == 2)
						isSuccess4 = true;
					else {
						System.out.println("\nZla wartosc - podaj wlasciwa: 1 - nie, 2 - tak\n");
					}
				}
				while(!isSuccess4);
				
				if(isRNA == 2) {
					
					seq1 = translateRNA(seq1, 1);
					seq2 = translateRNA(seq2, 2);
					System.out.println("Sekwencje RNA po translacji do aminokwasow:\n"+seq1+"\n"+seq2);
				}
				
		
		int [][] matrix = new int [seq2.length()][seq1.length()]; // rozszerzona macierz stworzona na podstawie macierzy z pliku aby latwiej liczylo sie dopasowania w metodach
		
		// rozszerzenie macierzy do dopasowania
		for(int i = 0; i < seq2.length(); i++)
		{		
		    for(int j = 0; j < seq1.length(); j++)
		    {
		    	
		    		    	
		    	
		        for(int k=0; k < rows; k++)
		        {
		        	for(int l=0; l < columns; l++)
		        	{
		        		if(seq1.charAt(j) == seq1Pom[l] && seq2.charAt(i) == seq2Pom[k])
		        		{
		        			matrix[i][j] = matrixRead[k][l];
		        		}
		        	}
		        }
		        

		        
		    }
		}
		
		//pobranie z macierzy kary za indel i wstawienie do Mapy
		Map<Character, Integer> gapPenalty = new HashMap<Character, Integer>(){};
		for(int i = 0; i < rows; i++)
		{		
		    for(int j = 0; j < columns; j++)
		    {
        		if(seq1Pom[j] == '-')
        		{
        			gapPenalty.put(seq2Pom[i],matrixRead[i][j]);
        		}
        		else if(seq2Pom[i] == '-')
        		{
        			gapPenalty.put(seq1Pom[j],matrixRead[i][j]);
        		}
        		
		    }
		}
		System.out.println("\nKara za indel: \n"+ gapPenalty);
		
		
		// wypisanie rozszerzonej macierzy
		System.out.println("\nMacierz rozszerzona: \n");
		System.out.print("\t");
				for(int i = 0; i < seq1.length(); ++i) {
						System.out.print(seq1.charAt(i)+"\t");
				}
				System.out.println("\n");
				
				for(int i = 0; i < seq2.length(); ++i)
				{
						System.out.print(seq2.charAt(i)+"\t");
					
				    for(int j = 0; j < seq1.length(); ++j)
				    {
				        
				    		System.out.print(matrix[i][j]+"\t");
				        
				    }
				    System.out.println();
				}
		
		//wprowadzanie typu macierzy z linii polecen
		boolean isSuccess2 = false;
		int matrixType = 0;
		System.out.println("\nWprowadz jaki typ macierzy uzyc: 1 - macierz odleg³osci, 2 - macierz podobienstwa\n");
		scan = null;
		do 
		{
			try {
				scan = new Scanner(System.in);
				matrixType = scan.nextInt();
			} catch (InputMismatchException e) {
				isSuccess2 = false;
			}
			
			if(matrixType == 1 || matrixType == 2)
				isSuccess2 = true;
			else {
				System.out.println("\nZla wartosc - podaj wlasciwa: 1 - macierz odleg³osci, 2 - macierz podobienstwa\n");
			}
		}
		while(!isSuccess2);
		
		//wprowadzanie typu dopasowania z linii polecen
		boolean isSuccess3 = false;
		int alignmentType = 0;
		System.out.println("\nWprowadz jaki typ dopasowania obliczyæ: 1 - globalne, 2 - lokalne\n");
		scan = null;
		do 
		{
			try {
				scan = new Scanner(System.in);
				alignmentType = scan.nextInt();
			} catch (InputMismatchException e) {
				isSuccess3 = false;
			}
			
			if(alignmentType == 1 || alignmentType == 2)
				isSuccess3 = true;
			else {
				System.out.println("\nZla wartosc - podaj wlasciwa: 1 - globalne, 2 - lokalne\n");
			}
		}
		while(!isSuccess3);
		
		
		SequenceAlignment sequenceAligner = new SequenceAlignment(gapPenalty);
		
		if(alignmentType == 1)
			sequenceAligner.calcGlobalAlignment(seq1, seq2, true, matrix, matrixType, isRNA);
		else if (alignmentType ==2)
			sequenceAligner.calcLocalAlignment(seq1, seq2, true, matrix, matrixType, isRNA);
		
		sc1.close();
		sc2.close();
		sc3.close();
		scan.close();
	}

}