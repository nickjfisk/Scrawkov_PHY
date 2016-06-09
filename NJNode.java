import java.util.ArrayList;

/**
 * A Node class to represent the NJ graph search.
 * A different node was needed than UPGMA as the 
 * boundary conditions and edge cases need be handled a little
 * differently.
 * @author J. Nick Fisk
 *
 */
public class NJNode {
	ArrayList<String> names=new ArrayList<String>();
	double dist1=Double.MIN_VALUE;
	double dist2=Double.MIN_VALUE;
	ArrayList<NJNode> njnodes=new ArrayList<NJNode>();
	String newick="";
	/**
	 * Constructor for the first nodes
	 * @param name
	 */
	public NJNode (String name){
		this.names.add(name);
		name=name.replaceAll("[()]", "");
		name=name.replaceAll(">", "");
		//name=name.replaceAll("\\","");
		//name=name.replaceAll("/","");
	//	if(name.length()>30){
		//	name=name.substring(0, 29);
		//}
		name=name.replaceAll(" ", "_");
		name=name.replaceAll("-", "_");
		name.replaceAll(",", "_");
		this.newick=name;
	}
	/**
	 * Constructor for grouping together two nodes 
	 * with two different lengths.
	 * @param n1
	 * @param n2
	 * @param dist1
	 * @param dist2
	 */
	public NJNode(NJNode n1, NJNode n2, double dist1, double dist2){
		this.njnodes.add(n1);
		this.njnodes.add(n2);
		this.names.addAll(n1.names);
		this.names.addAll(n2.names);
		this.newick="("+n1.newick+":"+dist1+","+n2.newick+":"+dist2+")";
	}
	/**
	 * String representation of node
	 */
	public String toString(){
		return(String.valueOf(this.names));
	}
	
	@Override
	public boolean equals(Object o){
		if(o==this){
			return true;
		}
		if(!(o instanceof Node)){
			return false;
		}
		Node o2=(Node)o;
		for(String n: this.names){
			if(o2.names.contains(n)==false){
				return false;
			}
		}
		for(String n: o2.names){
			if(this.names.contains(n)==false){
				return false;
			}
		}
		return true;
	}
	/**
 	* needed for the equals override
 	*/
	public int hashCode(){
		return this.names.hashCode();
	}
}