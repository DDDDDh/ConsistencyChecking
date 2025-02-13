package src.cn.edu.nju.moon.consistency.datastructure;

import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.Map;

import src.cn.edu.nju.moon.consistency.model.GlobalData;
import src.cn.edu.nju.moon.consistency.model.operation.ReadIncOperation;

/**
 * @description for each {@link ReadIncOperation}, it keeps a record
 * 	of the latest preceding WRITE {@link ReadIncOperation} for each variable.
 * 
 * @author hengxin
 * @date 2012-12-9
 * 
 * @see ReadIncOperation#lwMap
 */
public class LatestWriteMap
{
	/**
	 * <var, operation> pair, in which:
	 * 		String: variable
	 * 		{@link ReadIncOperation}: latest WRITE
	 * 
	 * @constraints {@link ReadIncOperation} must be WRITE and perform on var
	 */
	private Map<String, ReadIncOperation> latestWriteMap = null;

	/**
	 * @constructor initialize {@link #latestWriteMap}
	 */
	public LatestWriteMap()
	{
		this.latestWriteMap = new HashMap<String, ReadIncOperation>();
		
		for (String var : GlobalData.VARSET)
			this.latestWriteMap.put(var, null);
	}
	
	/**
	 * update this {@link #latestWriteMap} according to @param riop
	 * 
	 * @param riop dependent {@link ReadIncOperation} to update this {@link #latestWriteMap} 
	 */
	public void updateLatestWrite(ReadIncOperation riop)
	{
		ReadIncOperation riop_var_wriop = null;
		for (String var : GlobalData.VARSET)	// variable by variable
		{
			riop_var_wriop = riop.getLatestWriteMap().getLatestWrite(var);
			if (riop_var_wriop != null) {
//				System.out.println("not a null" + riop.toString());
				if (this.getLatestWrite(var) == null || this.getLatestWrite(var).getWid() < riop_var_wriop.getWid())
					this.latestWriteMap.put(var, riop_var_wriop);    // update to the most latest one
			}
//			else{
//				System.out.println("a null");
//			}
		}
		
		if (riop.isWriteOp() && ! riop.getWritetoOrder().isEmpty())	// there is some READ reads from @param riop
		{
			String var = riop.getVariable();
			if (this.getLatestWrite(var) == null)
				this.latestWriteMap.put(var, riop);
			else if (this.getLatestWrite(var).getWid() < riop.getWid())
				this.latestWriteMap.put(var, riop);
		}
	}
	
	/**
	 * get the latest WRITE {@link ReadIncOperation} on variable @param var
	 * @param var variable
	 * @return latest WRITE {@link ReadIncOperation} on @param var
	 * 
	 * @warning You SHOULD check whether the result is NULL.
	 * TODO: refactor to deal with the "NULL" return value.
	 */
	public ReadIncOperation getLatestWrite(String var)
	{
		return this.latestWriteMap.get(var);
	}
	
	/**
	 * compare two WRITE {@link ReadIncOperation}s and return the most latest one
	 * 
	 * @param wriop1 one WRITE {@link ReadIncOperation}
	 * @param wriop2 one WRITE {@link ReadIncOperation}
	 * @return the most latest WRITE
	 * 
	 * @constraints @param wriop1 and @param wriop2 must perform on the same variable
	 * @deprecated not necessary
	 */
	private ReadIncOperation maxWrite(ReadIncOperation wriop1, ReadIncOperation wriop2)
	{
		assertTrue("Only WRITEs on the same variable can be compared", wriop1.isWriteOp() && wriop2.isWriteOp());
		
		if (wriop1.getWid() < wriop2.getWid())
			return wriop2;
		return wriop1;
	}
}
