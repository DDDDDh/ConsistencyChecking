package src.cn.edu.nju.moon.consistency.datastructure;

import static org.junit.Assert.assertTrue;

import java.util.Set;

import src.cn.edu.nju.moon.consistency.model.operation.BasicOperation;
import src.cn.edu.nju.moon.consistency.model.operation.ReadIncOperation;
import src.cn.edu.nju.moon.consistency.model.process.ReadIncProcess;

/**
 * @description the leftmost READ {@link ReadIncOperation} 
 * 		reachable from this {@link ReadIncOperation}; 
 * 
 * @author hengxin
 * @date 2012-12-9
 * 
 * @see {@link ReadIncOperation}{@link #earlistRead}
 */
public class EarliestRead
{
	private int earlistRead = -1;
	
	/**
	 * initialize the earliest {@link ReadIncOperation} (by {@link #earlistRead})
	 * reachable to be @param rriop
	 * 
	 * @param rriop the earliest reachable {@link ReadIncOperation}
	 * 
	 * @constraints @param rriop must be READ {@link ReadIncOperation}
	 */
	public void initEarlistRead(ReadIncOperation rriop)
	{
		assertTrue("Using READ ReadIncOperation as initial value", rriop.isReadOp());
		
		this.earlistRead = rriop.getEarliestRead().earlistRead;
	}
	
	/**
	 * update the earlist read {@link #earlistRead} according to @param wriop
	 * 
	 * @param riop used to update {@link #earlistRead};
	 * the new earlist read is the smaller between two
	 * @return old {@link #earlistRead}
	 * 
	 * @modified hengxin on 2013-1-8
	 * @reason   add the assertion code: toEarliestRead != -1;
	 */
	public int updateEarliestRead(ReadIncOperation riop)
	{
		int oldEarlistRead = this.earlistRead;
		int toEarlistRead = riop.getEarliestRead().getEarlistReadInt();
		if (toEarlistRead != -1 && this.earlistRead > toEarlistRead)	// take the smaller one
			this.earlistRead = toEarlistRead;
		
		return oldEarlistRead;
	}
	
	/**
	 * update the earlist read {@link #earlistRead} according to a list of WRITE {@link ReadIncOperation}s
	 * 
	 * @param wriopSet set of {@link ReadIncOperation}s
	 * @return old {@link #earlistRead}
	 */
	public int updateEarliestRead(Set<BasicOperation> wriopSet)
	{
		int oldEarlistRead = this.earlistRead;
		
		for (BasicOperation wriop : wriopSet)	// take the smallest one
			this.updateEarliestRead((ReadIncOperation) wriop);

		return oldEarlistRead;
	}
	
	/**
	 * identify the (W,R) pair with which Rule (c) W'WR order can be applied
	 * 
	 * @param oldEarlistRead @param wriop {@link #earlistRead} of (W')
	 * @param wriop W' in W'WR order
	 * @param master_proc master {@link ReadIncProcess}
	 * 
	 * @return W in W'WR order
	 * 
	 * @warning you SHOULD check whether the returned value is null
	 */
	public ReadIncOperation identify_wrPair(int oldEarlistRead, ReadIncOperation wriop, ReadIncProcess master_proc)
	{
		assertTrue("W' in W'WR order must be WRITE", wriop.isWriteOp());
		
		String var = wriop.getVariable();
		int newEarlistRead = wriop.getEarliestRead().earlistRead;
		ReadIncOperation rriop = null;
		// identify R: the first READ with var(R) = var && r.index is in [oldEarlistRead, newEarlistRead)
		for (int index = newEarlistRead; index < oldEarlistRead; index++)
		{
			rriop = (ReadIncOperation) master_proc.getOperation(index);
			if (rriop.isReadOp() && rriop.getVariable().equals(var))
				return (ReadIncOperation) rriop.getReadfromWrite();	// return the dictating WRITE
		}
		
		return null;
	}
	
	public int getEarlistReadInt()
	{
		return this.earlistRead;
	}
	
	public void setEarlistReadInt(int er)
	{
		this.earlistRead = er;
	}
}
