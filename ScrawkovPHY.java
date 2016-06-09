import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
/**
 * 
 * @author J. Nick Fisk
 * Main class for the alignment free phylogeny algorithm, Scrawkov-Phy.
 */
public class ScrawkovPHY{
	//all the data from one file (as separated by newlines)
	public static ArrayList<String> allData=new ArrayList<String>();
	//the names of the taxa in order of appearance (to conserve sequence, name order)
	public static ArrayList<String> finchInOrder=new ArrayList<String>();
	//sequences in name-sequencce order in the desired format 
	public static ArrayList<String> formatedSeqList=new ArrayList<String>();
	//ID to bird object pairings. Basically ID sequence pairs. 
	public static HashMap<String, Bird> birdMap=new HashMap<String, Bird>();
	//contains all of the standard codons. 
	public static ArrayList<String> allCodons=new ArrayList<String>();
	//contains all of the standard dinucs
	public static ArrayList<String> allDinucs=new ArrayList<String>();
	//key for going from codons to AA
	public static HashMap<String,String> codonAAMap=new HashMap<String,String>();
	//goes from AA to a list of possible codons
	public static HashMap<String,ArrayList<String>> AAcodonMap=new HashMap<String, ArrayList<String>>();
	//has the scores for one gene tree.
	public static HashMap<String,HashMap<String,Double>>initScores=new HashMap<String,HashMap<String,Double>>();
	//public static HashMap<String,HashMap<Integer,HashMap<String,Double>>> QHMMTri=new HashMap<String,HashMap<Integer,HashMap<String,Double>>>();
	//public static HashMap<String,HashMap<Integer,HashMap<String,Double>>> QHMMDi=new HashMap<String,HashMap<Integer,HashMap<String,Double>>>();
	public static HashMap<String,HashMap<String, Double>> QHMMTri_final=new HashMap<String,HashMap<String, Double>>();
	public static HashMap<String,HashMap<String, Double>> QHMMDi_final=new HashMap<String,HashMap<String, Double>>();
	public static HashMap<String,HashMap<String,HashMap<String,Double>>> allScoresByGene= new HashMap<String,HashMap<String,HashMap<String,Double>>>();
	
	public static int k=-1; //kmer size. If negative (default), two will be used. 
	private static boolean doNJ=true; //NJ tree desired
	private static String njOut=""; //name of njOut file
	public static boolean maxByPair=true; //default behaviour should be true
	private static boolean normalize=true; //default should be true
	private static String geneOrSpecies="gene"; //default to gene
	private static ArrayList<String> fastaFiles=new ArrayList<String>(); //keep track of all the files for species trees
	private static String outputFile=null; //UPGMA output
	private static String usage="Usage: java ScrawkovPHY <treetype> <inputFasta> <optionalParams>";
	private static double weightGC=10; 
	private static double weightBias=.01;
	private static double weightTriFreq=1;
	private static double weightDiFreq=1;
	private static double weightTriQHMM=.00005;
	private static double weightDiQHMM=.00005;
	private static boolean outputAll=false;
	private static boolean bootSeq=false;
	private static String outBoot=null;
	private static boolean diag=false;
	

	//Globally available parameters
	/**
	 * Parses the command line arguments and sets globals accordingly
	 * @param params are the command line arguments in a arraylist
	 */
	private static void parseParams(ArrayList<String> params){
		//ask the magic counch
		//"Nothing."
		//THE COUNCH HAS SPOKEN!
		if(params.get(0).equals("gene")||params.get(0).equals("Gene")){
			geneOrSpecies="gene";
			params.remove(0);
			fastaFiles.add(params.get(0)); //only one fasta file
			params.remove(0);
			File f = new File(fastaFiles.get(0));
			if(!f.exists()){
				System.err.println("File does not exist");
				System.err.println(usage);
				System.exit(1);
			}
			if(f.isDirectory()){
				System.err.println("Input file is actually directory");
				System.err.println(usage);
				System.exit(1);
			}
		}
		else if(params.get(0).equals("species")||params.get(0).equals("Species")){
			geneOrSpecies="species";
			params.remove(0);
			if(params.isEmpty()){
				System.err.println("No input files provided");
				System.err.println("Exiting...");
				System.err.println(usage);
				System.exit(1);
			}
			if(params.get(0).equals("--folder")){
				File folder=new File(params.get(1));
				File[] listOfFiles=folder.listFiles();
				for(int i=0; i<listOfFiles.length;i++){
					if(listOfFiles[i].isFile()){
						fastaFiles.add(listOfFiles[i].getAbsolutePath());
					}
				}
				params.remove(0);
				params.remove(0);
			}
			String curFile=params.get(0);
			while((curFile.charAt(0)=='-')==false){
				fastaFiles.add(curFile);
				params.remove(0);
				curFile=params.get(0);
			}
			for(String i:fastaFiles){
				File f = new File(i);
				if(!f.exists()){
					System.err.println("File does not exist");
					System.err.println(usage);
					System.exit(1);
				}
				if(f.isDirectory()){
					System.err.println("Input file is actually directory");
					System.err.println(usage);
					System.exit(1);
				}
			}
			///merge the gene and species stuff
		}
		else if(params.get(0).equals("--help")||params.get(0).equals("-h")){
			printHelp();
			System.exit(0);
		}
		else{
			System.err.println("Invalid tree type. Valid values are 'gene' and 'species'");
			System.err.println(usage);
			System.exit(1);
		}
		while(!params.isEmpty()){
			String curFlag=params.get(0);
			params.remove(0);
			if(params.isEmpty()){
				System.err.println("Unpaired flag "+ curFlag);
				System.err.println("Exiting...");
				System.err.println(usage);
				System.exit(1);
			}
			String curValue=params.get(0);
			params.remove(0);
			switch(curFlag){
				case "--maxByPair":
					switch(curValue){
						case "T":
							maxByPair=true;
							break;
						case "t":
							maxByPair=true;
							break;
						case "true":
							maxByPair=true;
							break;
						case "True":
							maxByPair=true;
							break;
						case "TRUE":
							maxByPair=true;
							break;
						case "F":
							maxByPair=false;
							break;
						case "f":
							maxByPair=false;
							break;
						case "false":
							maxByPair=false;
							break;
						case "False":
							maxByPair=false;
							break;
						case "FALSE":
							maxByPair=false;
							break;
						default:
							System.err.println("Invalid value for maxByPair");
							System.err.println("Exiting...");
							System.err.println(usage);
							System.exit(1);
							break;
					}
					break;
				case "--outputFile":
					outputFile=curValue;
					break;
				case "-o":
					print(curValue);
					outputFile=curValue;
					break;
				case "--kSize":
					k=Integer.parseInt(curValue);
					break;
				case "--doNJ":
					switch(curValue){
						case "T":
							doNJ=true;
							break;
						case "t":
							doNJ=true;
							break;
						case "true":
							doNJ=true;
							break;
						case "True":
							doNJ=true;
							break;
						case "TRUE":
							doNJ=true;
							break;
						case "F":
							doNJ=false;
							break;
						case "f":
							doNJ=false;
							break;
						case "false":
							doNJ=false;
							break;
						case "False":
							doNJ=false;
							break;
						case "FALSE":
							doNJ=false;
							break;
						default:
							System.err.println("Invalid value for doNJ");
							System.err.println("Exiting...");
							System.err.println(usage);
							System.exit(1);
							break;
						}
						break;
				case "--diagMat":
					switch(curValue){
						case "T":
							diag=true;
							break;
						case "t":
							diag=true;
							break;
						case "true":
							diag=true;
							break;
						case "True":
							diag=true;
							break;
						case "TRUE":
							diag=true;
							break;
						case "F":
							diag=false;
							break;
						case "f":
							diag=false;
							break;
						case "false":
							diag=false;
							break;
						case "False":
							diag=false;
							break;
						case "FALSE":
							diag=false;
							break;
						default:
							System.err.println("Invalid value for diagMat");
							System.err.println("Exiting...");
							System.err.println(usage);
							System.exit(1);
							break;
						}
						break;
				case "--njOut":
					njOut=curValue;
					break;
				case "--bootSeq":
					switch(curValue){
						case "T":
							bootSeq=true;
							break;
						case "t":
							bootSeq=true;
							break;
						case "true":
							bootSeq=true;
							break;
						case "True":
							bootSeq=true;
							break;
						case "TRUE":
							bootSeq=true;
							break;
						case "F":
							bootSeq=false;
							break;
						case "f":
							bootSeq=false;
							break;
						case "false":
							bootSeq=false;
							break;
						case "False":
							bootSeq=false;
							break;
						case "FALSE":
							bootSeq=false;
							break;
						default:
							System.err.println("Invalid value for bootSeq");
							System.err.println("Exiting...");
							System.err.println(usage);
							System.exit(1);
							break;
						}
						break;	
					
				case "--outputAll":
					switch(curValue){
						case "T":
							outputAll=true;
							break;
						case "t":
							outputAll=true;
							break;
						case "true":
							outputAll=true;
							break;
						case "True":
							outputAll=true;
							break;
						case "TRUE":
							outputAll=true;
							break;
						case "F":
							outputAll=false;
							break;
						case "f":
							outputAll=false;
							break;
						case "false":
							outputAll=false;
							break;
						case "False":
							outputAll=false;
							break;
						case "FALSE":
							outputAll=false;
							break;
						default:
							System.err.println("Invalid value for outputAll");
							System.err.println("Exiting...");
							System.err.println(usage);
							System.exit(1);
							break;
						}
						break;
				case "--normalize":
					switch(curValue){
						case "T":
							normalize=true;
							break;
						case "t":
							normalize=true;
							break;
						case "true":
							normalize=true;
							break;
						case "True":
							normalize=true;
							break;
						case "TRUE":
							normalize=true;
							break;
						case "F":
							normalize=false;
							break;
						case "f":
							normalize=false;
							break;
						case "false":
							normalize=false;
							break;
						case "False":
							normalize=false;
							break;
						case "FALSE":
							normalize=false;
							break;
						default:
							System.err.println("Invalid value for normalize");
							System.err.println("Exiting...");
							System.err.println(usage);
							System.exit(1);
							break;
					}
					break;
					//--weightGC --weightCodon --weightTri --weightDi
				case "--weightGC":
					if(isDouble(curValue)){
						weightGC=Double.parseDouble(curValue);
					}
					else{
						System.err.println("Invalid value for weightGC");
						System.err.println("Exiting...");
						System.err.println(usage);
						System.exit(1);
					}
					break;
				case "--weightBias":
					if(isDouble(curValue)){
						weightBias=Double.parseDouble(curValue);
					}
					else{
						System.err.println("Invalid value for weightBias");
						System.err.println("Exiting...");
						System.err.println(usage);
						System.exit(1);
					}
					break;
				case "--weightTriQHMM":
					if(isDouble(curValue)){
						weightTriQHMM=Double.parseDouble(curValue);
					}
					else{
						System.err.println("Invalid value for weightTriQHMM");
						System.err.println("Exiting...");
						System.err.println(usage);
						System.exit(1);
					}
					break;
				case "--weightDiQHMM":
					if(isDouble(curValue)){
						weightDiQHMM=Double.parseDouble(curValue);
					}
					else{
						System.err.println("Invalid value for weightDiQHMM");
						System.err.println("Exiting...");
						System.err.println(usage);
						System.exit(1);
					}
					break;
				case "--weightTriFreq":
					if(isDouble(curValue)){
						weightTriFreq=Double.parseDouble(curValue);
					}
					else{
						System.err.println("Invalid value for weightTriFreq");
						System.err.println("Exiting...");
						System.err.println(usage);
						System.exit(1);
					}
					break;
				case "--weightDiFreq":
					if(isDouble(curValue)){
						weightDiFreq=Double.parseDouble(curValue);
					}
					else{
						System.err.println("Invalid value for weightDiFreq");
						System.err.println("Exiting...");
						System.err.println(usage);
						System.exit(1);
					}
					break;
				case "--help":
					printHelp();
					System.exit(0);
					break;
				default:
					System.err.println("Invalid flag");
					System.err.println("Exiting...");
					System.err.println(usage);
					System.exit(1);
					break;
			}
			
		}
	}

	/**
	 * The main function of the program, calling gene or species tree methods.
	 * It also instantiates all of the reference data structures. 
	 * @param args command line arguments
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		//populate some internal data structures
		popAllCodons();
		popAllDinucs();
		popCodonMap();
		popAAtoCodonMap();
		//get the parameters into an arrayList
		ArrayList<String> params=new ArrayList<String>();
		int parCount=0;
		while(parCount<args.length){
			params.add(args[parCount]);
			parCount+=1;
		}
		//parse the params
		parseParams(params);
		if(geneOrSpecies.equals("gene")||geneOrSpecies.equals("Gene")){
			doGeneTree(params,fastaFiles.get(0));//get fasta for now
			System.exit(0);
		}
		else if(geneOrSpecies.equals("species")||geneOrSpecies.equals("Species")){
			doSpeciesTree(params);
			ArrayList<String> allGeneTrees=new ArrayList<String>(allScoresByGene.keySet());
			ArrayList<String> allOrgs=new ArrayList<String>();
			for(String i: allGeneTrees){
				allOrgs.addAll(allScoresByGene.get(i).keySet());
			}
			Set<String> allOrgsSet=new HashSet<String>( allOrgs);
			allOrgs=new ArrayList<String>(allOrgsSet);
			HashMap<String,HashMap<String,Double>>finScores=new HashMap<String,HashMap<String,Double>>();
			HashMap<String,HashMap<String,Integer>> countMap=new HashMap<String,HashMap<String,Integer>>();
			for(String i: allOrgs){
				HashMap<String,Double>tempu=new HashMap<String,Double>();
				HashMap<String,Integer>tempu2=new HashMap<String,Integer>();
				for(String j: allOrgs){
					if(i.equals(j)==false){
						tempu.put(j, 0.0);
						tempu2.put(j,0);
					}
				}
				finScores.put(i, tempu);
				countMap.put(i,tempu2);
			}
			//only valid for species tree construction
			if(normalize==true){
				double max=getMaxScore();
				double miniMax=Double.MIN_VALUE;
				for(String i: allScoresByGene.keySet()){
					for(String j: allScoresByGene.get(i).keySet()){
						for(String k: allScoresByGene.get(i).get(j).keySet()){
							if(miniMax<allScoresByGene.get(i).get(j).get(k)){
								//find the largest element in the map
								miniMax=allScoresByGene.get(i).get(j).get(k);
							}
						}
					}
					//normalize everything to be proportional
					for(String j: allScoresByGene.get(i).keySet()){
						for(String k: allScoresByGene.get(i).get(j).keySet()){
							allScoresByGene.get(i).get(j).put(k, (allScoresByGene.get(i).get(j).get(k)*(max/miniMax)));
						}
					}
				}
			}
			//sum all the scores
			for(String j: allScoresByGene.keySet()){
				for(String q: allScoresByGene.get(j).keySet()){
					for(String h: allScoresByGene.get(j).get(q).keySet()){
						finScores.get(q).put(h, finScores.get(q).get(h)+allScoresByGene.get(j).get(q).get(h));
						countMap.get(q).put(h, countMap.get(q).get(h)+1);
					}
				}
			}
			//normalize by number of times a species appears.
			ArrayList<String>allSpeciesOccurances=new ArrayList<String>();
			for(String i: allScoresByGene.keySet()){	
				allSpeciesOccurances.addAll(allScoresByGene.get(i).keySet());
			}
			for(String i: countMap.keySet()){
				for(String j: countMap.get(i).keySet()){
					if(countMap.get(i).get(j)!=0){
						finScores.get(i).put(j, (finScores.get(i).get(j)/countMap.get(i).get(j)));
					}
				}
			}
			initScores=finScores;
			//build the UPGMA Tree
			HashMap<String,HashMap<String,Double>> wkMap=new HashMap<String,HashMap<String,Double>>(initScores);
			HashMap<Node,HashMap<Node,Double>>nodeDist=new HashMap<Node,HashMap<Node,Double>>();
			for(String s: initScores.keySet()){
				Node sn=new Node(s);
				nodeDist.put(sn, new HashMap<Node,Double>());
				for(String r: initScores.keySet()){
					if(s.equals(r)==false){
						Node rn=new Node(r);
						nodeDist.get(sn).put(rn, initScores.get(s).get(r));
					}
				}
			}
			ArrayList<String> allNames=new ArrayList<String>(wkMap.keySet());
			
			Node curMinNode1=null;
			Node curMinNode2=null;
			double curMinScore=Double.MAX_VALUE;
			HashMap<String,Boolean>inNode=new HashMap<String,Boolean>();
			for(String i: initScores.keySet()){
				inNode.put(i, false);
			}
			ArrayList<Node>nodeList=new ArrayList<Node>();
			while(nodeDist.keySet().size()>1){
				curMinScore=Double.MAX_VALUE;
				for(Node i: nodeDist.keySet()){
					for(Node ent:nodeDist.get(i).keySet()){
						if(nodeDist.get(i).get(ent)<curMinScore){
							curMinNode1=i;
							curMinNode2=ent;
							curMinScore=nodeDist.get(i).get(ent);
						}
					}
				}
				Node nw2=new Node((Node)curMinNode1,(Node)curMinNode2, curMinScore);
				nodeDist.remove(curMinNode1);
				nodeDist.remove(curMinNode2);
				for(Node n: nodeDist.keySet()){
					if(nodeDist.get(n).containsKey(curMinNode1)){
						nodeDist.get(n).remove(curMinNode1);
					}
					if(nodeDist.get(n).containsKey(curMinNode2)){
						nodeDist.get(n).remove(curMinNode2);
					}
				}
				nodeDist.put(nw2, new HashMap<Node,Double>());
				for(Node n: nodeDist.keySet()){
					if(n.equals(nw2)==false){
						//calc dist and put
						double calcDist=Node.calcDistance(nw2, n);
						nodeDist.get(nw2).put(n, calcDist);
						nodeDist.get(n).put(nw2, calcDist);
					}	
				}
			}
			for(Node fin: nodeDist.keySet()){
				print("Printing UPGMA Tree");
				print(fin.newick+";");
				if(outputFile!=null){
					PrintWriter out = new PrintWriter(outputFile+".newick");
					out.print(fin.newick+";");
					out.close();
				}
			}	
		}
	}
	
	/**
	 * Constructs a single gene tree from a single fasta file
	 * @param params, the command line arguments as an arraylist
	 * @param file, the fasta file
	 * @throws IOException
	 */
	public static void doGeneTree(ArrayList<String> params, String file) throws IOException{
		//initialize data structures
		ArrayList<String> allData=new ArrayList<String>();
		formatedSeqList=new ArrayList<String>();
		HashMap<String, Bird> birdMap=new HashMap<String, Bird>();
		initScores=new HashMap<String,HashMap<String,Double>>();;
		FileReader fr = new FileReader(new File(file)); 
		BufferedReader br = new BufferedReader(fr);  
		String line;
		//read in sequence
		while((line = br.readLine()) != null)
		{ 
		    line = line.trim(); // remove leading and trailing whitespace
		    allData.add(line);
		}
		String holdMe = "";
		for (String element : allData) {
			if (element.charAt(0) == '>') {
				if (holdMe != "") {
					formatedSeqList.add(holdMe);
					holdMe = "";
				}
				formatedSeqList.add(element);
			}
			else {
				holdMe += element;
			}
		}
		formatedSeqList.add(holdMe);
		fr.close();
		//make Bird Objects
		for(int i=0; i<formatedSeqList.size(); i++){
			String id=formatedSeqList.get(i);
			String seq=formatedSeqList.get(i+1);
			Bird newbie=new Bird(id,seq);
			birdMap.put(id,newbie);
			i+=1;
		}
		if(bootSeq==true){
			genBootstrapSeqs(true, 3, .85, birdMap);
		}
		//populate hidden markov model measures
		QHMMALL(birdMap);
		HashMap<String,ArrayList<String>>keepTrack=new HashMap<String,ArrayList<String>>();
		//score all
		for(String id: birdMap.keySet()){
			for(String id2:birdMap.keySet()){
				if(id.equals(id2)==false){
					if(keepTrack.containsKey(id)){
						if(keepTrack.get(id).contains(id2)==false){
							Double d=computeANDcombine(birdMap.get(id),birdMap.get(id2));
							if(initScores.containsKey(id)==false){
								initScores.put(id, new HashMap<String,Double>());
							}
							initScores.get(id).put(id2, d);
						}
					}
					else{
						keepTrack.put(id, new ArrayList<String>());
						keepTrack.get(id).add(id2);
						initScores.put(id, new HashMap<String, Double>());
						if(keepTrack.containsKey(id2)==false){
							keepTrack.put(id2, new ArrayList<String>());
							keepTrack.get(id2).add(id);
							Double d=computeANDcombine(birdMap.get(id),birdMap.get(id2));
							initScores.get(id).put(id2, d);
						}
					}
				}
			}
		}
		for(String i: initScores.keySet()){
			for(String j: initScores.keySet()){
				if(i.equals(j)==false){
					if(initScores.get(i).get(j)==null){
						initScores.get(i).put(j, initScores.get(j).get(i));
					}
				}
			}
		}
		/////CLUSTERING AND TREE BUILDING NEXT!/////
		HashMap<String,HashMap<String,Double>> wkMap=new HashMap<String,HashMap<String,Double>>(initScores);
		HashMap<Node,HashMap<Node,Double>>nodeDist=new HashMap<Node,HashMap<Node,Double>>();
		for(String s: initScores.keySet()){
			Node sn=new Node(s);
			nodeDist.put(sn, new HashMap<Node,Double>());
			for(String r: initScores.keySet()){
				if(s.equals(r)==false){
					Node rn=new Node(r);
					nodeDist.get(sn).put(rn, initScores.get(s).get(r));
				}
			}
		}
		toMatrix(nodeDist);
		Node curMinNode1=null;
		Node curMinNode2=null;
		double curMinScore=Double.MAX_VALUE;
		while(nodeDist.keySet().size()>1){
			curMinScore=Double.MAX_VALUE;
			for(Node i: nodeDist.keySet()){
				for(Node ent:nodeDist.get(i).keySet()){
					if(nodeDist.get(i).get(ent)<curMinScore){
						curMinNode1=i;
						curMinNode2=ent;
						curMinScore=nodeDist.get(i).get(ent);
					}
				}
			}
			Node nw2=new Node((Node)curMinNode1,(Node)curMinNode2, curMinScore);
			nodeDist.remove(curMinNode1);
			nodeDist.remove(curMinNode2);
			for(Node n: nodeDist.keySet()){
				if(nodeDist.get(n).containsKey(curMinNode1)){
					nodeDist.get(n).remove(curMinNode1);
				}
				if(nodeDist.get(n).containsKey(curMinNode2)){
					nodeDist.get(n).remove(curMinNode2);
				}
			}
			nodeDist.put(nw2, new HashMap<Node,Double>());
			for(Node n: nodeDist.keySet()){
				if(n.equals(nw2)==false){
					//calc dist and put
					double calcDist=Node.calcDistance(nw2, n);
					nodeDist.get(nw2).put(n, calcDist);
					nodeDist.get(n).put(nw2, calcDist);
				}
			}
		}
		for(Node fin: nodeDist.keySet()){
			print("Printing UPGMA Tree");
			print(fin.newick+";");
			if(outputFile!=null){
				if(geneOrSpecies.equals("gene")){
					PrintWriter out = new PrintWriter(outputFile);
					out.print(fin.newick+";");
					out.close();
				}
			}
			if(outputAll==true){
				if(geneOrSpecies.equals("species")){
					File f = new File("geneTreesForSpeciesTree");
					if (f.exists() && f.isDirectory()) {
					   PrintWriter out=new PrintWriter(file+"gene_tree.newick");
					   out.print(fin.newick);
					   out.close();
					}
					else{
						f.mkdir();
						PrintWriter out=new PrintWriter(f+file+"_gene_tree.newick");
						out.print(fin.newick);
						out.close();
					}
				}
			}
		}
	}
	/**
	 * Formats a HashMap of Distances to a PHYLIP formatted distance matrix
	 * @param nodeDist
	 */
	private static void toMatrix(HashMap<Node, HashMap<Node, Double>> nodeDist) {		
		ArrayList<Node> nodes=new ArrayList<Node>(nodeDist.keySet());
		int matSize=nodes.size();
		double [][] scoreMatrix=new double[matSize][matSize];
		ArrayList<String> stringsToWrite=new ArrayList<String>();
		ArrayList<String> namesList=new ArrayList<String>();
		for(Node n : nodes){
			String temp=n.names.get(0).substring(1, n.names.get(0).length());
			namesList.add(temp);
			temp=temp+"          ";
			if(temp.length()>10){
				temp=temp.substring(0,10);
			}
			stringsToWrite.add(temp+" ");
		}
		try {
			PrintWriter write=new PrintWriter("distanceMatrix.txt","UTF-8");
			write.println(" "+matSize);
			//i is row and j is col
			for(int i=0; i< nodes.size();i++){
				Node n=nodes.get(i);
				for(int j=0; j<nodes.size();j++){
					if(i==j){
						scoreMatrix[i][j]=0.0;
					}
					else{
						Node n2=nodes.get(j);
						scoreMatrix[i][j]=Node.calcDistance(n, n2);
					}
				}
			}
			//diagonal formmat of phylip distance matrix
			if(diag){
				int count=matSize;
				for(int i=0; i<matSize;i++){
					String tmpStr=stringsToWrite.get(i);
					for(int j=0+count;j<matSize;j++){
						if(j!=matSize-1){
							tmpStr=tmpStr+scoreMatrix[i][j-count]+" ";
						}
						else{
							tmpStr=tmpStr+scoreMatrix[i][j-count];
						}
					}
					count=count-1;
					write.println(tmpStr);
				}
			}
			else{
				for(int i=0;i<matSize;i++){
					String tempStr=stringsToWrite.get(i);
					for(int j=0;j<matSize;j++){
						if(j!=matSize-1){
							tempStr=tempStr+scoreMatrix[i][j]+" ";
						}
						else{
							write.println(tempStr+scoreMatrix[i][j]);
						}
					}
				}
			}
			write.close();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(doNJ==true){
			doNJTree(scoreMatrix,nodeDist,matSize,namesList);
		}
	}
/**
 * Constructs a neighbor joining tree in newick format
 * @param scoreMatrix
 * @param nodeDist
 * @param matSize
 * @param namesList
 */
	public static void doNJTree(double[][] scoreMatrix, HashMap<Node, HashMap<Node, Double>> nodeDist,int matSize,ArrayList<String> namesList){
		//i is row, j is column
		HashMap<NJNode,HashMap<NJNode,Double>> NJNodeMap=new HashMap<NJNode,HashMap<NJNode,Double>>();
		ArrayList<NJNode> NJList=new ArrayList<NJNode>();
		ArrayList<ArrayList<Double>> M=new ArrayList<ArrayList<Double>>();
		ArrayList<ArrayList<Double>> Mp=new ArrayList<ArrayList<Double>>();
		for(int i=0; i<namesList.size(); i++){
			NJNode tmp=new NJNode(namesList.get(i));
			NJList.add(i,tmp);
		}
		for(int i=0; i<NJList.size();i++){
			for(int j=0;j<NJList.size();j++){
				if(i!=j){
					if(NJNodeMap.containsKey(NJList.get(i))==false){
						NJNodeMap.put(NJList.get(i), new HashMap<NJNode,Double>());
					}
					if(NJNodeMap.get(NJList.get(i)).containsKey(NJList.get(j))==false){
						NJNodeMap.get(NJList.get(i)).put(NJList.get(j), scoreMatrix[i][j]);
					}
				}
			}
		}
		while(NJNodeMap.keySet().size()>2){
			HashMap<NJNode,HashMap<NJNode,Double>> Q=new HashMap<NJNode,HashMap<NJNode,Double>>();
			HashMap<NJNode,Double> T=new HashMap<NJNode,Double>();
			//compute T, the sum total of rows
			double totalRow=0;
			for(NJNode i: NJNodeMap.keySet()){
				for(NJNode j: NJNodeMap.keySet()){
					if(!i.equals(j)){
						totalRow+=NJNodeMap.get(i).get(j);
					}
				}
				T.put(i, totalRow);
				totalRow=0;
			}
			NJNode curMinNode1=null;
			NJNode curMinNode2=null;
			double curMinScore=Double.POSITIVE_INFINITY;
			//compute Q, the transformed matrix for neighbor joining minimal selection
			for(NJNode i: NJNodeMap.keySet()){
				Q.put(i, new HashMap<NJNode,Double>());
				for(NJNode j:NJNodeMap.keySet()){
					if(!i.equals(j)){
						if(Q.containsKey(j)==false){
							Q.put(j, new HashMap<NJNode, Double>());
						}
						double tmp=(((NJNodeMap.keySet().size()-2)*NJNodeMap.get(i).get(j))-T.get(i)-T.get(j));
						Q.get(i).put(j,tmp);
						Q.get(j).put(i, tmp);
						if(tmp<curMinScore){
							if(!i.equals(j)){
								curMinNode1=i;
								curMinNode2=j;
								curMinScore=tmp;
							}
						}
					}
				}
			}
			Q.clear();
			//Find branch lengths
			double deltaIJ=((T.get(curMinNode1)-T.get(curMinNode2))/(NJNodeMap.size()-2));
			double LLi=0.5*(NJNodeMap.get(curMinNode1).get(curMinNode2)+deltaIJ);
			double LLj=0.5*(NJNodeMap.get(curMinNode1).get(curMinNode2)-deltaIJ);
			//make new node
			NJNode tempuNode=new NJNode(curMinNode1,curMinNode2,LLi,LLj);
			//compute new distances to other nodes and add to map.
			double joinDist=0;
			HashMap<NJNode, Double> joinMap=new HashMap<NJNode,Double>();
			for(NJNode n: NJNodeMap.keySet()){
				if(!n.equals(curMinNode1)){
					if(!n.equals(curMinNode2)){
						joinDist=(NJNodeMap.get(n).get(curMinNode1)+NJNodeMap.get(n).get(curMinNode2)-NJNodeMap.get(curMinNode1).get(curMinNode2))/2;
						joinMap.put(n, joinDist);
					}
				}
			}
			//remove the two old nodes...
			NJNodeMap.remove(curMinNode1);
			NJNodeMap.remove(curMinNode2);
			for(NJNode n: NJNodeMap.keySet()){
				if(NJNodeMap.get(n).containsKey(curMinNode1)){
					NJNodeMap.get(n).remove(curMinNode1);
				}
				if(NJNodeMap.get(n).containsKey(curMinNode2)){
					NJNodeMap.get(n).remove(curMinNode2);
				}
			}
			//add new nodes
			for(NJNode n: NJNodeMap.keySet()){
				NJNodeMap.get(n).put(tempuNode, joinMap.get(n));
			}
			NJNodeMap.put(tempuNode, joinMap);
		}
		ArrayList<NJNode>tmp2=new ArrayList<NJNode>(NJNodeMap.keySet());
		NJNode fin=tmp2.get(0);
		NJNode fin2=tmp2.get(1);
		print("Printing NJ Tree");
		String newickNJTree="("+fin2.newick+":"+NJNodeMap.get(fin).get(fin2)+","+fin.newick+");";
		print(newickNJTree);
		if(njOut.equals("")==false){
			PrintWriter write;
			try {
				write = new PrintWriter(njOut,"UTF-8");
				write.write(newickNJTree);
				write.close();
			} catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
			}
		}
		else{
			PrintWriter write;
			try {
				write = new PrintWriter("NJ_Tree.nwk","UTF-8");
				write.write(newickNJTree);
				write.close();
			} catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
			}			
		}
	}

	//START SPECIES TREE CONSTRUCTION (resumed in main for adaptive scope reasons)////
	/**
	 * Scores every gene tree needed for the construction of a species tree
	 * by calling the doGeneTree method and storing the result in a global hashmap
	 * @param params parameters in ArrayList that will affect the settings of the tree building algorithm
	 * @throws IOException
	 */
	public static void doSpeciesTree(ArrayList<String>params) throws IOException{
		for(String name : fastaFiles){
			doGeneTree(params,name);
			allScoresByGene.put(name, initScores);
		}
		for(String o : allScoresByGene.keySet()){
			HashMap<String, HashMap<String,Double>>temp=allScoresByGene.get(o);
			print(temp.keySet());
			print("");
		}
	}
	/**
	 * Finds the largest score of any tree for any gene
	 * Used only in species tree method with scaling enabled
	 * @return a double of the highest score observed
	 */
	public static double getMaxScore(){
		double curMax=Double.MIN_VALUE;
		for(String i: allScoresByGene.keySet()){
			for(String j: allScoresByGene.get(i).keySet()){
				for(String k: allScoresByGene.get(i).get(j).keySet()){
					if(allScoresByGene.get(i).get(j).get(k)>curMax){
						curMax=allScoresByGene.get(i).get(j).get(k);
					}
				}
			}
		}
		return curMax;
	}
	
	/**
	 * Calculates the scores for each of the components of
	 * the algorithm and combines them to get a composite score
	 * @param bird1 taxa 1 to get the score between
	 * @param bird2 taxa 2 to get the score between
	 * @return the composite score
	 */
	public static double computeANDcombine(Bird bird1, Bird bird2){
		return ((scoreBirdCodon(bird1,bird2)*weightTriFreq)+(scoreBirdDinuc(bird1,bird2)*weightDiFreq)+
		+(scoreGC(bird1,bird2)*weightGC)+(scoreCodonBias(bird1,bird2)*weightBias)+
		(QHMMTri_final.get(bird1.id).get(bird2.id)*weightTriQHMM)+(QHMMDi_final.get(bird1.id).get(bird2.id)*weightDiQHMM));		
	}
	
	/**
	 * Computes the score for the difference in codon usage (codon bias)
	 * between two taxa.
	 * @param bird1 taxa 1
	 * @param bird2 taxa 2
	 * @return the score, stored as a double
	 */
	public static double scoreCodonBias(Bird bird1, Bird bird2){
		double totalScore=0.0;
		double protScore=0.0;
		int numPotCodons=0;
		double scoreLogCheck=0.0;
		for(String aa1:AAcodonMap.keySet()){
			protScore=0;
			if(bird1.biasMap.containsKey(aa1)==false){
				//bird one doesnt, but bird 2 does
				if(bird2.biasMap.containsKey(aa1)==true){
					numPotCodons=bird2.biasMap.get(aa1).keySet().size();
					for(String cod:bird2.biasMap.get(aa1).keySet()){
						protScore+=Math.abs(log2(bird2.biasMap.get(aa1).get(cod)/numPotCodons));
					}
					totalScore+=protScore;
				}
				//neither do
				else{
					//pass, rather, add 0 to total.
				}
			}
			//bird1 has the aa1
			else{
				//bird 2 doesn't
				if(bird2.biasMap.containsKey(aa1)==false){
					numPotCodons=bird1.biasMap.get(aa1).keySet().size();
					for(String cod:bird1.biasMap.get(aa1).keySet()){
						protScore+=Math.abs(log2(bird1.biasMap.get(aa1).get(cod)/numPotCodons));
					}
					totalScore+=protScore;
				}
				//they both have aa of interest
				else{
					Set<String> unionCodons=new HashSet<String>();
					unionCodons.addAll(bird1.biasMap.get(aa1).keySet());
					unionCodons.addAll(bird2.biasMap.get(aa1).keySet());
					numPotCodons=unionCodons.size();
					for(String cod:unionCodons){
						//bird 2 doesn't have the codon
						if(bird2.biasMap.get(aa1).containsKey(cod)==false){
							//that means bird one must!
							protScore+=Math.abs(log2(bird1.biasMap.get(aa1).get(cod)/numPotCodons));
						}
						//bird two does have the codon
						else{
							//check if bird one does...
							//no?
							if(bird1.biasMap.get(aa1).containsKey(cod)==false){
								protScore+=Math.abs(log2(bird2.biasMap.get(aa1).get(cod)/numPotCodons));
							}
							//yes, they both do
							else{
								scoreLogCheck=Math.abs(bird1.biasMap.get(aa1).get(cod)-bird2.biasMap.get(aa1).get(cod)/numPotCodons);
								if(scoreLogCheck!=0){
									protScore+=Math.abs(log2(scoreLogCheck));
								}
								else{
									//add zero
								}
							}
						}
					}
					totalScore+=protScore;
				}
			}
		}
		totalScore=totalScore/20;
		return totalScore;
	}

	/**
	 * Calculates the score for a difference in codon transition probabilities between
	 * two taxa for a given gene.
	 * @param bird1 taxa 1
	 * @param bird2 taxa 2
	 * @return score, stored as a double
	 */
	public static double scoreBirdCodon(Bird bird1, Bird bird2){
		double bestDist=Double.MAX_VALUE;
		for(int i: bird1.codonFrequencies.keySet()){
			double score1=0.0;
			double score2=0.0;
			double distScore=0.0;
			for(String c1:allCodons){
				for(String c2:allCodons){
					if(bird1.codonFrequencies.get(i).containsKey(c1)==false){
						score1=0.0;
					}
					else{
						if(bird1.codonFrequencies.get(i).get(c1).containsKey(c2)==false){
							score1=0.0;
						}
						else{
							score1=bird1.codonFrequencies.get(i).get(c1).get(c2);
						}
					}
					if(bird2.codonFrequencies.get(i).containsKey(c1)==false){
						score2=0.0;
					}
					else{
						if(bird2.codonFrequencies.get(i).get(c1).containsKey(c2)==false){
							score2=0.0;
						}
						else{
							score2=bird2.codonFrequencies.get(i).get(c1).get(c2);
						}
					}
					distScore+=Math.abs(score1-score2);
				}
			}
			if(distScore<bestDist){
				bestDist=distScore;
			}
		}
		return bestDist;
	}
	
	/**
	 * Constructs and computes the scores for the QHMM for both Tri and Di nucleotides
	 * The results are not returned--rather, stored as a feature in the inputted birdMap.
	 * @param birdMap, the internal data structure used to calculate QHMM
	 */
	public static void QHMMALL(HashMap<String, Bird> birdMap){
		HashMap<String,HashMap<Integer,HashMap<String,Double>>> QHMMTri=new HashMap<String,HashMap<Integer,HashMap<String,Double>>>();
		HashMap<String,HashMap<Integer,HashMap<String,Double>>> QHMMDi=new HashMap<String,HashMap<Integer,HashMap<String,Double>>>();
		int[] iterlist=new int[3];
		iterlist[0]=1;
		iterlist[1]=2;
		iterlist[2]=3;
		for(String name: birdMap.keySet()){
			HashMap<Integer,HashMap<String,Double>> tempu=new HashMap<Integer,HashMap<String,Double>>();
			for(int i: iterlist){	
				HashMap<String,Double>tempu2=new HashMap<String,Double>();
					for(String name2: birdMap.keySet()){
						if(name.equals(name2)==false){
							tempu2.put(name2, 0.0);
						}
					tempu.put(i, tempu2);
					}
				}
			QHMMTri.put(name, tempu);
		}
		for(String name: QHMMTri.keySet()){
			Bird tempuBird=birdMap.get(name);
			String tseq=tempuBird.sequence;
			for(int i: iterlist){
				while(tseq.length()>5){
					String codon1=tseq.substring(0, 3);
					String codon2=tseq.substring(3, 6);
					for(String name2: birdMap.keySet()){
						if(!name2.equals(name)){
							if(birdMap.get(name2).codonFrequencies.get(i).containsKey(codon1)){
								if(birdMap.get(name2).codonFrequencies.get(i).get(codon1).containsKey(codon2)){
									QHMMTri.get(name).get(i).put(name2, (QHMMTri.get(name).get(i).get(name2))+(Math.abs(log2(birdMap.get(name2).codonFrequencies.get(i).get(codon1).get(codon2)))));
								}	
							}
						}
					}
					tseq=tseq.substring(3);
				}
				tseq=tempuBird.sequence;
			}
			for(String name2: QHMMTri.keySet()){
				if(name.equals(name2)==false){
					double score1=QHMMTri.get(name).get(1).get(name2);
					double score2=QHMMTri.get(name).get(2).get(name2);
					double score3=QHMMTri.get(name).get(3).get(name2);
					double bestScore=Math.min(Math.min(score1, score2), score3);
					if(QHMMTri_final.containsKey(name)==false){
						HashMap<String,Double> ntempu=new HashMap<String,Double>();
						ntempu.put(name2, bestScore);
						QHMMTri_final.put(name, ntempu);
					}
					else{
						QHMMTri_final.get(name).put(name2, bestScore);
					}
				}
			}
		}		
		int numDo=2;
		if(k>1){
			numDo=k;
		}
		iterlist=new int[numDo];
		for(int m=0;m<numDo;m++){
			iterlist[m]=m+1;
		}
		for(String name: birdMap.keySet()){
			HashMap<Integer,HashMap<String,Double>> tempu=new HashMap<Integer,HashMap<String,Double>>();
			for(int i: iterlist){	
				HashMap<String,Double>tempu2=new HashMap<String,Double>();
					for(String name2: birdMap.keySet()){
						if(name.equals(name2)==false){
							tempu2.put(name2, 0.0);
						}
					tempu.put(i, tempu2);
					}
				}
			QHMMDi.put(name, tempu);
		}
		for(String name: QHMMDi.keySet()){
			Bird tempuBird=birdMap.get(name);
			String tseq=tempuBird.sequence;
			int whileLen=(2*numDo)-1;
			for(int i: iterlist){
				while(tseq.length()>whileLen){
					String codon1=tseq.substring(0, numDo);
					String codon2=tseq.substring(numDo, numDo+numDo);
					for(String name2: birdMap.keySet()){
						if(!name2.equals(name)){
							if(birdMap.get(name2).dinucFrequencies.get(i).containsKey(codon1)){
								if(birdMap.get(name2).dinucFrequencies.get(i).get(codon1).containsKey(codon2)){
									QHMMDi.get(name).get(i).put(name2, (QHMMDi.get(name).get(i).get(name2))+(Math.abs(log2(birdMap.get(name2).dinucFrequencies.get(i).get(codon1).get(codon2)))));
								}	
							}
						}
					}
					tseq=tseq.substring(3);
				}
				tseq=tempuBird.sequence;
			}
			for(String name2: QHMMDi.keySet()){
				if(name.equals(name2)==false){
					double score1=QHMMDi.get(name).get(1).get(name2);
					double score2=QHMMDi.get(name).get(2).get(name2);
					double bestScore=Math.min(score1, score2);
					if(QHMMDi_final.containsKey(name)==false){
						HashMap<String,Double> ntempu=new HashMap<String,Double>();
						ntempu.put(name2, bestScore);
						QHMMDi_final.put(name, ntempu);
					}
					else{
						QHMMDi_final.get(name).put(name2, bestScore);
					}
				}
			}
		}
	}
	
	/**
	 * Calculates the score for a difference in dinucleotide transition probabilities
	 * between two taxa for a given gene.
	 * @param bird1 taxa 1
	 * @param bird2 taxa 2
	 * @return score, stored as a double
	 */
	public static double scoreBirdDinuc(Bird bird1, Bird bird2) {
		double score1 = 0.0;
		double score2 = 0.0;
		double distScore = 0.0;
		for (int i : bird1.dinucFrequencies.keySet()) {
			for (String c1 : allDinucs) {
				for (String c2 : allDinucs) {
					if (bird1.dinucFrequencies.containsKey(c1) == false) {
						score1 = 0.0;
					} else {
						if (bird1.dinucFrequencies.get(c1).containsKey(c2) == false) {
							score1 = 0.0;
						} else {
							score1 = bird1.dinucFrequencies.get(i).get(c1)
									.get(c2);
						}
					}
					if (bird2.dinucFrequencies.containsKey(c1) == false) {
						score2 = 0.0;
					} else {
						if (bird2.dinucFrequencies.get(i).get(c1)
								.containsKey(c2) == false) {
							score2 = 0.0;
						} else {
							score2 = bird2.dinucFrequencies.get(i).get(c1)
									.get(c2);
						}
					}
					distScore += Math.abs(score1 - score2);
				}
			}
		}
		return distScore;
	}
	/**
	 * A shortcut function that gets the log
	 * base two on the double passed in
	 * @param num number to get the log base 2 of
	 * @return the log base two of the number
	 */
	
	public static double log2(double num){
		return Math.log(num)/Math.log(2);
	}
	
	/**
	 * Returns the difference in GC scores between two taxa
	 * @param bird1 taxa 1
	 * @param bird2 taxa 2
	 * @return score, stored as a double
	 */
	public static double scoreGC(Bird bird1, Bird bird2){
		return Math.abs(bird1.GC_score-bird2.GC_score);
	}
	
	/**
	 * A shortcut function for System.out.println
	 * @param o, an object to print
	 */
	public static void print(Object o){
		System.out.println(o);
	}
	/**
	 * Make a file
	 * @param filenameANDpath, the absolute path of the desired new file
	 */
	public static void mkFile(File filenameANDpath){
		File newfile=filenameANDpath;
		try{
			newfile.createNewFile();
		}
		catch(IOException ioe){
			System.err.println("Error in making files\n "+ioe);
		}
	}
	
	/**
	 * Makes subsequence files for use in the bootstrapping method.
	 * @param sameStart, a boolean indicating if all sequences should be sub-setted starting at the same location
	 * @param numBoots, the number of times to bootstrap
	 * @param propBases, the proportion of bases wanted in the bootstrap
	 * @param birdMap, the data structure containing all the sequences
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 */
	public static void genBootstrapSeqs(boolean sameStart, int numBoots, double propBases, HashMap<String, Bird> birdMap) throws FileNotFoundException, UnsupportedEncodingException{
		new File("tmp").mkdir();
		new File("bootTrees").mkdir();
		int numBases=0;
		int startIndex=0;
		int sizeOfSeq=0;
		for(int i=0; i<numBoots; i++){
			if(sameStart==false){
				try {
					PrintWriter writer=new PrintWriter("tmp/bootstrapIter"+i+".fasta","UTF-8");
					for(String berd: birdMap.keySet()){
						String id=birdMap.get(berd).id;
						String seq=birdMap.get(berd).sequence;
						writer.println(id);
						sizeOfSeq=seq.length();
						numBases=(int) Math.round(propBases*sizeOfSeq);
						startIndex=genRandInt(0,(sizeOfSeq-numBases)-2);
						writer.println(seq.substring(startIndex, startIndex+numBases-1));
						print(seq.substring(startIndex, startIndex+numBases-1));
					}
					writer.close();
				}
			catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
				}
			}
			else{
				String brd=(String) birdMap.keySet().toArray()[0];
				String id=birdMap.get(brd).id;
				String seq=birdMap.get(brd).sequence;
				sizeOfSeq=seq.length();
				numBases=(int) Math.round(propBases*sizeOfSeq);
				startIndex=genRandInt(0,(sizeOfSeq-numBases)-2);
				PrintWriter writer=new PrintWriter("tmp/bootstrapIter"+i+".fasta","UTF-8");
				for(String berd: birdMap.keySet()){
					id=birdMap.get(berd).id;
					seq=birdMap.get(berd).sequence;
					writer.println(id);
					if(startIndex+numBases>seq.length()){
						writer.println(seq.substring(startIndex, seq.length()-1));
					}
					else{
						writer.println(seq.substring(startIndex, startIndex+numBases-1));
					}
					
				}
				writer.close();
			}
		}
	}
	/**
	 * Shortcut function for generating a random integer.
	 * @param Min the mimimum bound
	 * @param Max the maximum bound
	 * @return a randomly generated integer
	 */
	public static int genRandInt(int Min, int Max){
		return (Min + (int)(Math.random() * ((Max - Min) + 1)));
	}
	
	/**
	 * Just see if a value is parse-able to a double
	 * Used in command line argument processing. 
	 * @param value the string being checked
	 * @return a boolean representing saying if it is a double
	 */
	public static boolean isDouble(String value) {
	    try {
	        Double.parseDouble(value);
	        return true;
	    } catch (NumberFormatException e) {
	        return false;
	    }
	}
	
	/**
	 * Internally used command that populates
	 * the codon to amino acid hashmap
	 */
	public static void popCodonMap(){
		codonAAMap.put("TTT","F");
		codonAAMap.put("TTC","F");
		codonAAMap.put("TTA","L");
		codonAAMap.put("TTG","L");
		codonAAMap.put("CTT","L");
		codonAAMap.put("CTC","L");
		codonAAMap.put("CTG","L");
		codonAAMap.put("CTA","L");
		codonAAMap.put("ATT","I");
		codonAAMap.put("ATC","I");
		codonAAMap.put("ATA","I");
		codonAAMap.put("ATG","M");
		codonAAMap.put("GTT","V");
		codonAAMap.put("GTC","V");
		codonAAMap.put("GTA","V");
		codonAAMap.put("GTG","V");
		codonAAMap.put("TCT","S");
		codonAAMap.put("TCC","S");
		codonAAMap.put("TCA","S");
		codonAAMap.put("TCG","S");
		codonAAMap.put("CCT","P");
		codonAAMap.put("CCC","P");
		codonAAMap.put("CCA","P");
		codonAAMap.put("CCG","P");
		codonAAMap.put("ACT","T");
		codonAAMap.put("ACC","T");
		codonAAMap.put("ACA","T");
		codonAAMap.put("ACG","T");
		codonAAMap.put("GCT","A");
		codonAAMap.put("GCC","A");
		codonAAMap.put("GCA","A");
		codonAAMap.put("GCG","A");
		codonAAMap.put("TAT","Y");
		codonAAMap.put("TAC","Y");
		codonAAMap.put("TAA","STOP");
		codonAAMap.put("TAG","STOP");
		codonAAMap.put("CAT","H");
		codonAAMap.put("CAC","H");
		codonAAMap.put("CAA","Q");
		codonAAMap.put("CAG","Q");
		codonAAMap.put("AAT","N");
		codonAAMap.put("AAC","N");	
		codonAAMap.put("AAA","K");
		codonAAMap.put("AAG","K");
		codonAAMap.put("GAT","D");
		codonAAMap.put("GAC","D");		
		codonAAMap.put("GAA","E");
		codonAAMap.put("GAG","E");
		codonAAMap.put("TGT","C");
		codonAAMap.put("TGC","C");
		codonAAMap.put("TGA","STOP");
		codonAAMap.put("TGG","W");
		codonAAMap.put("CGT","R");
		codonAAMap.put("CGC","R");
		codonAAMap.put("CGA","R");
		codonAAMap.put("CGG","R");
		codonAAMap.put("AGT","S");
		codonAAMap.put("AGC","S");
		codonAAMap.put("AGA","R");
		codonAAMap.put("AGG","R");
		codonAAMap.put("GGT","G");
		codonAAMap.put("GGC","G");
		codonAAMap.put("GGA","G");
		codonAAMap.put("GGG","G");
	}
	/**
	 * Internally used function that builds the Amino acid
	 * to Codon map used in algorithm. It relies on popCodonMap()
	 * to have been called prior to the calling of this method.
	 */
	public static void popAAtoCodonMap(){
		for(String i:codonAAMap.keySet()){
			String AA=codonAAMap.get(i);
			if(AAcodonMap.containsKey(AA)==false){
				ArrayList<String>tempu=new ArrayList<String>();
				tempu.add(i);
				AAcodonMap.put(AA,tempu);
			}
		}
	}
	/**
	 * A shortcut method for getting translating codons into amino acids
	 * @param codon the codon to be translated
	 * @return the 1 letter amino acid code as String
	 */
	public static String translateCodon(String codon){
		return codonAAMap.get(codon);
	}
	/**
	 * Internally used command that populates a list of
	 * every possible dinucleotide combination
	 */
	public static void popAllDinucs(){
		allDinucs.add("AA");
		allDinucs.add("AT");
		allDinucs.add("AG");
		allDinucs.add("AC");
		allDinucs.add("TT");
		allDinucs.add("TA");
		allDinucs.add("TC");
		allDinucs.add("TG");
		allDinucs.add("GG");
		allDinucs.add("GC");
		allDinucs.add("GA");
		allDinucs.add("GT");
		allDinucs.add("CC");
		allDinucs.add("CA");
		allDinucs.add("CT");
		allDinucs.add("CG");
	}
	/**
	 * Internally used command that populates a list of every possible
	 * codon.
	 */
	public static void popAllCodons(){
		allCodons.add("TTT");
		allCodons.add("TTC");
		allCodons.add("TTG");
		allCodons.add("TTA");
		allCodons.add("TCT");
		allCodons.add("TCC");
		allCodons.add("TCA");
		allCodons.add("TCG");
		allCodons.add("TAT");
		allCodons.add("TAC");
		allCodons.add("TAA");
		allCodons.add("TAG");
		allCodons.add("TGT");
		allCodons.add("TGC");
		allCodons.add("TGA");
		allCodons.add("TGG");
		allCodons.add("CTT");
		allCodons.add("CTC");
		allCodons.add("CTA");
		allCodons.add("CTG");
		allCodons.add("CCT");
		allCodons.add("CCC");
		allCodons.add("CCA");
		allCodons.add("CCG");
		allCodons.add("CAT");
		allCodons.add("CAC");
		allCodons.add("CAA");
		allCodons.add("CAG");
		allCodons.add("CGT");
		allCodons.add("CGC");
		allCodons.add("CGA");
		allCodons.add("CGG");
		allCodons.add("ATT");
		allCodons.add("ATC");
		allCodons.add("ATA");
		allCodons.add("ATG");
		allCodons.add("ACT");
		allCodons.add("ACC");
		allCodons.add("ACA");
		allCodons.add("ACG");
		allCodons.add("AAT");
		allCodons.add("AAC");
		allCodons.add("AAA");
		allCodons.add("AAG");
		allCodons.add("AGT");
		allCodons.add("AGC");
		allCodons.add("AGA");
		allCodons.add("AGG");
		allCodons.add("GTT");
		allCodons.add("GTC");
		allCodons.add("GTA");
		allCodons.add("GTG");
		allCodons.add("GCT");
		allCodons.add("GCC");
		allCodons.add("GCA");
		allCodons.add("GCG");
		allCodons.add("GAT");
		allCodons.add("GAC");
		allCodons.add("GAA");
		allCodons.add("GAG");
		allCodons.add("GGT");
		allCodons.add("GGC");
		allCodons.add("GGA");
		allCodons.add("GGG");
	}
	/**
	 * prints the help statement seen with --help or -h
	 */
	public static void printHelp(){
		print("---ScrawkovPHY.java---");
		print(usage);
		print("Example: java ScrawkovPHY gene exampleFasta.fasta --normalize T");
		print("Example: java ScrawkovPHY species exampleFastaForGene1.fasta exampleFastaForGene2.fasta --maxByPair T");
		print(" ");
		print("-------Flags available, with desciption of fucntion---------");
		print(" ");
		print("--folder");
		print("Use a folder rather than a series of files in species tree contsruction");
		print(" ");
		print("--doNJ");
		print("   Construct a neighbor joining tree in addition to the default UPGMA tree");
		print("   Default behaviour is true");
		print(" ");
		print("--njOut");
		print("   The file name for the NJ tree output. Only used if doNJ is true");
		print("   The default behaviour is to use the name njOut.nwk");
		print(" ");
		print("--maxByPair:");
		print("   maxByPair takes either true or false as a value");
		print("   If true, the algorithm will find and use the optimal reading frame for each pairwise sequence");
		print("   The default behaviour is true ");
		print(" ");
		print("--normalize:");
		print("   normalize takes either true or false as a value");
		print("   This parameter only affects the construction of species trees");
		print("   If true, the results of each gene tree are scaled to the range of the most disparate");
		print("   gene tree. That is, the one with the biggest difference between highest and lowest score");
		print("   The default behaviour is true ");
		print(" ");
		print("--outputFile");
		print("   outputFile takes the name of the file you want the tree to be output in");
		print("   The default behaviour is to only print the output to the terminal, not to a file");
		print("   It has the shortcut -o");
		print(" ");
		print("--kSize");
		print("   override the kmer size for dinucleotide aspect of QHMM.");
		print("   default value is -1, which is ignored by program.");
		print(" ");
		print("--bootSeq");
		print("   Experimental: Should subsequences of DNA be obtained for bootstrapping? Must be run with gene");
		print("   Default behaviour is false");
		print(" ");
		print("--diagMat");
		print("   Should the PHYLIP formatted distance matrix be written along the diagonal?");
		print("   If false, the full matrix will be written.");
		print("   The default behaviour is true.");
		print(" ");
		print("--bootOut");
		print("   Experimental: Name of outPut bootstrapped tree folder. If a value is entered, will use the species tree infrastructure to make");
		print("   bootstrap trees from the output of bootSeq. Must be run with species");
		print("   The default behaviour is to not make these trees.");
		print(" ");
		print("--help");
		print("   Print this help statement");
		print(" ");
		print("--weightGC, --weightBias, --weightTriFreq, --weightDiFreq, --weightTriQHMM, --weightDiQHMM");
		print("   These flags accept numeric input (decimals are fine)");
		print("   These arguments are the weights applied to the 6 features used to create the MCCI, which is a surrogate distance.");
		print("   They affect the GC content, codon bias, trinucleotide frequency, dinucleotide frequency,");
		print("    trinucleotide transistion QHMM, and dinucleotide transtition QHMM respectively.");
		print("   The default values are 10,0.1,1,1,0.00005,and 0.00005 respectively.");
	}
}
