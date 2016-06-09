import java.util.HashMap;
/**
 * Bird is a generic representation of a taxon and its sequence.
 * @author J. Nick Fisk
 *
 */
public class Bird {
	public double GC_score; //GC content proportion
	public double CodonFreqScore; 
	public String sequence; //store the sequence in RAM
	public String id; //name of sequence
	public HashMap<Integer,HashMap<String, HashMap<String, Double>>> codonFrequencies; //initMaps
	public HashMap<Integer, HashMap<String, HashMap<String, Double>>> dinucFrequencies; 
	public HashMap<String,HashMap<String,Double>> biasMap;
	
	/**
	 * An object representing the a taxon in the algorithm.
	 * It is names bird as it was developed on birds originally. 
	 * @param id the id of the taxon from the fasta file
	 * @param sequence the nucleotide sequence
	 */
	public Bird(String id, String sequence){
		this.GC_score=0.0;
		this.CodonFreqScore=0.0;
		this.sequence=sequence;
		this.id=id;
		this.dinucFrequencies=new HashMap<Integer,HashMap<String,HashMap<String, Double>>> ();
		this.codonFrequencies=new HashMap<Integer,HashMap<String,HashMap<String, Double>>> ();
		this.codonFrequencies=newCalcCodonFreq();
		//check k from main to see if we are doing generic length kmer rather than a hard 2.
		if(ScrawkovPHY.k<1){
			this.dinucFrequencies=newCalcDiFreq();
		}
		else{
			this.dinucFrequencies=calcKFreq(ScrawkovPHY.k);
		}
		this.GC_score=calcGC();
		this.biasMap=mkBiasMap();
	}
	/**
	 * Measure the features ultimately used in the calculation of codon bias
	 * @return HashMap of measured features
	 */
	public HashMap<String, HashMap<String,Double>> mkBiasMap(){
		String inter=this.sequence; //shallow copy of the sequence
		HashMap<String,HashMap<String,Double>> biasMap=new HashMap<String,HashMap<String,Double>>();
		String pIter=""; //holder variable
		//go through and chop down the sequence, translating each codon at the end. 
		while(inter.length()>3){
				pIter=ScrawkovPHY.translateCodon(inter.substring(0,3));
				if(biasMap.containsKey(pIter)==false){
					HashMap<String,Double> bcodonMap=new HashMap<String,Double>();
					bcodonMap.put(inter.substring(0,3), (double) 1);
					biasMap.put(pIter, bcodonMap);
				}
				else{
					if(biasMap.get(pIter).containsKey(inter)==false){
						biasMap.get(pIter).put(inter.substring(0,3), 1.0);
					}
					else{
						biasMap.get(pIter).put(inter.substring(0,3),biasMap.get(pIter).get(inter)+1);
					}
				}
				inter=inter.substring(3);
		}
		for(String aa: biasMap.keySet()){
			double count=0;
			for(String entry: biasMap.get(aa).keySet()){
				count+=biasMap.get(aa).get(entry);
			}
			for(String entry:biasMap.get(aa).keySet()){
				biasMap.get(aa).put(entry,biasMap.get(aa).get(entry)/count);
			}
		}
		return biasMap;
	}
	/**
	 * Measures trinucleotide occurrence frequency, as well as the transition frequency 
	 * in sequences.
	 * @param seq sequence being observed
	 * @return HashMap of measured features
	 */
	public HashMap<String, HashMap<String, Double>> calcCodonHelper(String seq){
		String tseq=seq;
		HashMap<String, HashMap<String, Double>> tmap=new HashMap<String, HashMap<String, Double>>();
		int numCodons=0;
		//crawl down and chop the sequence into smaller pieces of size 3.
		while(tseq.length()>5){
			numCodons+=1;
			String codon1=tseq.substring(0, 3);
			String codon2=tseq.substring(3, 6);
			if(!tmap.containsKey(codon1)){
				HashMap<String,Double>tempu=new HashMap<String, Double>();
				Double count=1.0;
				tempu.put(codon2,  count);
				tmap.put(codon1, tempu);
			}
			else{
				if(!tmap.get(codon1).containsKey(codon2)){
					tmap.get(codon1).put(codon2, (double)1);
				}
				else{
					tmap.get(codon1).put(codon2, tmap.get(codon1).get(codon2)+1);
				}
			}
			tseq=tseq.substring(3);
		}
		for(String map: tmap.keySet()){
			for(String entry: tmap.get(map).keySet()){
				tmap.get(map).put(entry, tmap.get(map).get(entry)/numCodons);
			}
		}
		return tmap;
				
	}
	/**
	 * Calls the helper function a differning number of times depending on the status 
	 * of the maxByPair flag. 
	 * @return HashMap of measured features
	 */
	public HashMap<Integer,HashMap<String, HashMap<String, Double>>> newCalcCodonFreq(){
		String tseq=this.sequence;
		HashMap<Integer,HashMap<String, HashMap<String, Double>>> tmap=new HashMap<Integer,HashMap<String, HashMap<String, Double>>>();
		int numToDo;
		//if maxByPair is true, each starting frame will be represented.
		if(ScrawkovPHY.maxByPair==true){
			numToDo=3;
		}
		else{
			numToDo=1;
		}
		int count=0;
		while(count<=numToDo){
			tmap.put(count, calcCodonHelper(tseq));
			tseq=tseq.substring(1,tseq.length());
			count+=1;
		}
		return tmap;
	}
	/**
	 * Measures dinucleotide frequencies and transition frequencies in a sequence
	 * @param seq the sequence being observed
	 * @return HashMap of measured features
	 */
	public HashMap<String, HashMap<String, Double>> calcDiHelper(String seq){
		String tseq=seq;
		HashMap<String, HashMap<String, Double>> tmap=new HashMap<String, HashMap<String, Double>>();
		int numCodons=0;
		while(tseq.length()>5){
			numCodons+=1;
			String codon1=tseq.substring(0, 2);
			String codon2=tseq.substring(2, 4);
			if(!tmap.containsKey(codon1)){
				HashMap<String,Double>tempu=new HashMap<String, Double>();
				Double count=1.0;
				tempu.put(codon2,  count);
				tmap.put(codon1, tempu);
			}
			else{
				if(!tmap.get(codon1).containsKey(codon2)){
					tmap.get(codon1).put(codon2, (double)1);
				}
				else{
					tmap.get(codon1).put(codon2, tmap.get(codon1).get(codon2)+1);
				}
			}
			tseq=tseq.substring(2);
		}
		for(String map: tmap.keySet()){
			for(String entry: tmap.get(map).keySet()){
				tmap.get(map).put(entry, tmap.get(map).get(entry)/numCodons);
			}
		}
		return tmap;

	}
	/**
	 * Depending on the value of maxByPair flag, calls the helper function a variable number of times.
	 * @return A HashMap of observed features
	 */
	public HashMap<Integer,HashMap<String, HashMap<String, Double>>> newCalcDiFreq(){
		String tseq=this.sequence;
		HashMap<Integer,HashMap<String, HashMap<String, Double>>> tmap=new HashMap<Integer,HashMap<String, HashMap<String, Double>>>();
		int numToDo;
		if(ScrawkovPHY.maxByPair==true){
			numToDo=2;
		}
		else{
			numToDo=1;
		}
		int count=0;
		while(count<=numToDo){
			tmap.put(count, calcDiHelper(tseq));
			tseq=tseq.substring(1,tseq.length());
			count+=1;
		}
		return tmap;
	}
	/**
	 * Performs the same function as calcCodonHelper and calcDinucHelper, but with a generic sized kmer.
	 * @param seq The sequence
	 * @param k the size of the kmer to be used
	 * @return
	 */
	public HashMap<String, HashMap<String, Double>> calcKHelper(String seq, int k){
		String tseq=seq;
		HashMap<String, HashMap<String, Double>> tmap=new HashMap<String, HashMap<String, Double>>();
		int numCodons=0;
		int whileLen= (k*2)+1;
		while(tseq.length()>whileLen){
			numCodons+=1;
			String codon1=tseq.substring(0, k);
			String codon2=tseq.substring(k, k+k);
			if(!tmap.containsKey(codon1)){
				HashMap<String,Double>tempu=new HashMap<String, Double>();
				Double count=1.0;
				tempu.put(codon2,  count);
				tmap.put(codon1, tempu);
			}
			else{
				if(!tmap.get(codon1).containsKey(codon2)){
					tmap.get(codon1).put(codon2, (double)1);
				}
				else{
					tmap.get(codon1).put(codon2, tmap.get(codon1).get(codon2)+1);
				}
			}
			tseq=tseq.substring(k);
		}
		for(String map: tmap.keySet()){
			for(String entry: tmap.get(map).keySet()){
				tmap.get(map).put(entry, tmap.get(map).get(entry)/numCodons);
			}
		}
		return tmap;

	}
	/**
	 * Depending on the value of maxByPair flag, calls the helper function a variable number of times.
	 * @return A HashMap of observed features
	 */
	public HashMap<Integer,HashMap<String, HashMap<String, Double>>> calcKFreq(int k){
		String tseq=this.sequence;
		HashMap<Integer,HashMap<String, HashMap<String, Double>>> tmap=new HashMap<Integer,HashMap<String, HashMap<String, Double>>>();
		int numToDo;
		if(ScrawkovPHY.maxByPair==true){
			numToDo=k;
		}
		else{
			numToDo=1;
		}
		int count=0;
		while(count<=numToDo){
			tmap.put(count, calcKHelper(tseq,ScrawkovPHY.k));
			tseq=tseq.substring(1,tseq.length());
			count+=1;
		}
		return tmap;
	}
	/**
	 * Measures the GC content of a sequence. 
	 * @return the difference in gc content as a double
	 */
	public double calcGC(){
		String tseq=this.sequence;
		double count = tseq.length() - tseq.replace("G", "").length();
		count+=tseq.length()-tseq.replace("C", "").length();
		return count/tseq.length();
	}
}
