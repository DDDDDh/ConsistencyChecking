package src.cn.edu.nju.moon.consistency.model.observation;

import java.util.Set;

import src.cn.edu.nju.moon.consistency.checker.ReadIncChecker;
import src.cn.edu.nju.moon.consistency.datastructure.GlobalActiveWritesMap;
import src.cn.edu.nju.moon.consistency.model.process.ReadIncProcess;
import src.cn.edu.nju.moon.consistency.ui.DotUI;

/**
 * @description filtered RawObservation:
 * 	all the READ operations are ignored except those reside on process with masterPid 
 * 
 * @author hengxin
 * @date 2012-12-7
 *
 */
public class ReadIncObservation extends BasicObservation
{
	/**
	 *  pool of all WRITE operations
	 *  in the PRAM-distinct problem, WRITE operations are distinct.
	 *  
	 * @modified hengxin on 2013-1-9
	 * @modification remove public static field {@link #WRITEPOOL}; 
	 * 	change it into protected; extracted to super class {@link BasicObservation}
	 */
//	public static Map<String, ReadIncOperation> WRITEPOOL = null;	
	
//	// {@link ReadIncProcess} with masterPid is to be checked against PRAM consistency
//	private int masterPid = -1;
	
	// global active WRITEs for each variable; used by {@link ReadIncChecker}
	private GlobalActiveWritesMap globalActiveWritesMap = new GlobalActiveWritesMap();
	
	/**
	 * get the filtered observation in which READ operations
	 * are ignored except those reside on process with @param masterPid
	 * 
	 * @param masterPid process with masterPid is kept the same
	 * @param rob RawObservation to be filtered
	 * 
	 * @see {@link ReadIncProcess}
	 */
	public ReadIncObservation(int masterPid, BasicObservation rob)
	{
		/** initialize to avoid {@link NullPointerException} */
		super(rob.getProcNum());
		
		this.masterPid = masterPid;
		
		for (int pid : rob.getProcMap().keySet())
			this.procMap.put(pid, new ReadIncProcess(masterPid, rob.getProcess(pid)));
		
		this.storeWrite2Pool();
		
		// ui
		DotUI.getInstance().visual_ob(this);
	}

	/**
	 * @return the {@link ReadIncProcess} with masterPid
	 * 	which is also the process to be checked against PRAM consistency
	 * 
	 * @see {@link ReadIncChecker}
	 */
	public ReadIncProcess getMasterProcess()
	{
		return (ReadIncProcess) this.procMap.get(this.masterPid);
	}
	
	/**
	 * @return {@link #globalActiveWritesMap}
	 */
	public GlobalActiveWritesMap getGlobalActiveWritesMap()
	{
		return this.globalActiveWritesMap;
	}
	
//	/**
//	 * preprocessing the {@link ReadIncObservation}, including
//	 * (1) establishing "program order" between {@link ReadIncOperation}
//	 * (2) establishing "write to order" between {@link ReadIncOperation}
//	 */
//	public void preprocessing()
//	{
//		this.establishProgramOrder();
//		this.establishWritetoOrder();
//	}
	
//	/**
//	 * does some READ {@link ReadIncOperation} read value from later WRITE
//	 * {@link ReadIncOperation} on the same {@link ReadIncProcess} with masterPid
//	 * 
//	 * @return true, if some READ {@link ReadIncOperation} read value from later WRITE
//	 * 	{@link ReadIncOperation} on the same {@link ReadIncProcess} with masterPid;
//	 * 		   false, o.w..
//	 * 
//	 * @see ReadIncProcess#readLaterWrite()
//	 * @see ReadIncChecker#check_part()
//	 */
//	public boolean readLaterWrite()
//	{
//		return ((ReadIncProcess) this.procMap.get(this.masterPid)).readLaterWrite();
//	}
//	
//	/**
//	 * is there no {@link ReadIncOperation} at all in the {@link ReadIncProcess}
//	 * with masterPid. if it is the case, PRAM consistency is satisfied trivially.
//	 * 
//	 * @return true, if no {@link ReadIncOperation} 
//	 * 	in {@link ReadIncProcess} with masterPid;	false, o.w..
//	 */
//	public boolean nullCheck()
//	{
//		if (this.procMap.get(this.masterPid).size() == 0)
//			return true;
//		return false;
//	}
	
//	/**
//	 * establishing "program order" between {@link ReadIncOperation}s
//	 */
//	private void establishProgramOrder()
//	{
//		Set<Integer> pids = this.procMap.keySet();
//		for (int pid : pids)
//			((ReadIncProcess) this.procMap.get(pid)).establishProgramOrder();
//	}
	
//	/**
//	 * establishing "write to order" between {@link ReadIncOperation}s 
//	 * and set rid for READ {@link ReadIncOperation} and wid for corresponding
//	 * WRITE {@link ReadIncOperation} 
//	 */
//	private void establishWritetoOrder()
//	{
//		// FIXME: NullPointer exception
//		
//		// all READ {@link ReadIncOperation}s are in the {@link ReadIncProcess} with #masterPid
//		((ReadIncProcess) this.procMap.get(this.masterPid)).establishWritetoOrder();
//	}
}
