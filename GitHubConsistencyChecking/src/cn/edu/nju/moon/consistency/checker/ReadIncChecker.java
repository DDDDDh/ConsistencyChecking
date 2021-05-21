package src.cn.edu.nju.moon.consistency.checker;

import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

//import com.sun.tools.javac.code.Lint;
import src.cn.edu.nju.moon.consistency.datastructure.GlobalActiveWritesMap;
import src.cn.edu.nju.moon.consistency.model.GlobalData;
import src.cn.edu.nju.moon.consistency.model.observation.BasicObservation;
import src.cn.edu.nju.moon.consistency.model.observation.ReadIncObservation;
import src.cn.edu.nju.moon.consistency.model.operation.BasicOperation;
import src.cn.edu.nju.moon.consistency.model.operation.RawOperation;
import src.cn.edu.nju.moon.consistency.model.operation.ReadIncOperation;
import src.cn.edu.nju.moon.consistency.model.process.ReadIncProcess;
import src.cn.edu.nju.moon.consistency.schedule.ISchedule;
import src.cn.edu.nju.moon.consistency.ui.DotUI;
//import sun.awt.image.ImageWatched;

/**
 * @author hengxin
 * 
 * @modified hengxin on 2013-1-5
 * @reason 	 using Template Method design pattern
 * @see 	 Checker
 */
public class ReadIncChecker extends Checker
{
	private ReadIncObservation riob = null;	/** {@link ReadIncObservation} with respect to some process to check **/
	
	/**
	 * Constructor
	 * @param riob {@link BasicObservation} to check
	 */
	public ReadIncChecker(BasicObservation rob)
	{
		super(rob);
	}
	
	/**
	 * Constructor
	 * @param riob	{@link BasicObservation} to check
	 * @param name	for {@link DotUI}; the name of file for visualization
	 */
	public ReadIncChecker(BasicObservation rob, String name)
	{
		super(rob, name);
	}
	
	/**
	 * Constructor
	 * @param riob	{@link BasicObservation} to check
	 * @param name	for {@link DotUI}; the name of file for visualization
	 * @param s 	record of checking results
	 */
	public ReadIncChecker(BasicObservation rob, String name, ISchedule s)
	{
		super(rob, name, s);
	}
	
	/**
	 * @return {link ReadIncObservation} with respect to @param masterPid to check
	 */
	@Override
	protected BasicObservation getMasterObservation(int masterPid)
	{
		return new ReadIncObservation(masterPid, super.rob);
	}
	
	/**
	 * check whether {@link #riob} satisfies PRAM Consistency:
	 * 
	 * @return true; if {@link #riob} satisfies PRAM Consistency; false, otherwise.
	 */
	@Override
	protected boolean check_part(BasicObservation rob)
	{
		assertTrue("check ReadIncObservation", rob instanceof ReadIncObservation);
		this.riob = (ReadIncObservation) rob;
		
		ReadIncProcess master_proc = this.riob.getMasterProcess();
		int master_size = master_proc.getOpNum();
		ReadIncOperation master_cur_rriop = null;
		BasicOperation bop = null;	
		boolean consistent = true;
		
		for (int index = 0; index < master_size; index++)
		{
			bop = master_proc.getOperation(index);
//			System.out.println("bop:" + bop.toString());
			if (bop.isReadOp())	// "ReadIncremental" checking algorithm is READ centric.
			{
//                System.out.println("rop:" + bop.toString());
                master_cur_rriop = (ReadIncOperation) bop;
//				System.out.println("bop earlist read:" + ((ReadIncOperation) bop).getEarliestRead().getEarlistReadInt());
				master_proc.set_cur_rriop(master_cur_rriop);	// set the current {@link ReadIncOperation} to check
				
				// (1) compute global active WRITEs
				boolean dominated = this.compute_globalActiveWritesMap(master_proc, master_proc.get_pre_rriop(), master_cur_rriop);
				// ui: for GAWM after computing intervals
				DotUI.getInstance().addGAWM(this.riob.getGlobalActiveWritesMap(), master_cur_rriop.toString(), 1);
				
				// (2) master_cur_rriop must read value from dw; apply Rule (c): W'WR order
				/**
				 * @modified hengxin on 2013-1-9
				 * @reason   if no other WRITEs than dw itself, it is not necessary to reschedule at all
				 */
				ReadIncOperation dw_cur_rriop = (ReadIncOperation) master_cur_rriop.getReadfromWrite();
//				if (this.reschedule_assertion(master_cur_rriop, dw_cur_rriop))
//				{
					if (this.readFromDW(master_cur_rriop, dw_cur_rriop))
					{
						consistent = false;	/** cycle **/
						System.out.println("cycle detected at riop:" + master_cur_rriop.toStringCom());
						break;
					}
					
					// (3) reschedule operations in r'-downset (i.e., master_pre_rriop-downset)
					if (dominated) {
                        this.loopTime++;
//                        System.out.println("loop inc for r':" + master_proc.get_pre_rriop() + " r:" + master_cur_rriop);
                        if (this.reschedule(master_cur_rriop)) {
							System.out.println("cycle detected at riop:" + master_cur_rriop.toStringCom() +" after reschedule");
                            consistent = false;    /** cycle **/
                            break;
                        }
                    }
//                    else{
//                        System.out.println("loop stay for r':" + master_proc.get_pre_rriop() + " r:" + master_cur_rriop);
//                    }
//				}
//				else{
//				    this.skipLoop++;
//                    System.out.println("do not need to reschedule by cur_rriop:" + master_cur_rriop + " dw:" + dw_cur_rriop + " process" + master_proc.getPid());
//                }
//				 ui: for GAWM after rescheduling (even if no rescheduling at all: for test)
				DotUI.getInstance().addGAWM(this.riob.getGlobalActiveWritesMap(), master_cur_rriop.toString(), 2);
				
				master_proc.advance_pre_rriop(master_cur_rriop);	// iterate over the next (R,R) pair
			}
		}
		
		// ui
		DotUI.getInstance().execute("readinc/" + name + "_" + master_proc.getPid());

        System.out.println("Loop time for inc alg of process " + riob.getMasterPid() + ": " + this.loopTime);

        return consistent;	/** no cycle; satisfying PRAM Consistency **/
	}

	/**
	 * computation of global active WRITEs map: 
	 * 		{@link GlobalActiveWritesMap} in {@link ReadIncObservation}
	 *  
	 * @param master_proc {@link ReadIncProcess} to be checked; it is the process with master_pid
	 * @param master_pre_rriop previous READ {@link ReadIncOperation} 
	 * @param master_cur_rriop current READ {@link ReadIncOperation} in consideration
	 * @return true, if "dw" is in master_pre_rriop-downset; false, otherwise.
	 * 
	 * @constraints both @param master_pre_rriop and @param master_cur_rriop must be READ {@link ReadIncOperation}
	 * 
	 * NOTE: the returned value indicates whether "dw" is in master_pre_rriop-downset;
	 * 	it is by-product of this method
	 */
	private boolean compute_globalActiveWritesMap(ReadIncProcess master_proc, ReadIncOperation master_pre_rriop, ReadIncOperation master_cur_rriop)
	{
		assertTrue("READ incremental: two arguments should be READ", master_pre_rriop.isReadOp() && master_cur_rriop.isReadOp());

//		System.out.println("Dealing with master_proc:" + master_proc.getPid() + " pre r:" + master_pre_rriop.toString() + " cur r:"  + master_cur_rriop.toString());
		
		boolean dominated = false;
		
		int master_pre_rindex = master_pre_rriop.getIndex();
		int master_cur_rindex = master_cur_rriop.getIndex();
		
		// (1) dealing with rr_interval: (master_pre_rriop, master_cur_rriop)
		ReadIncOperation pre_riop = master_pre_rriop;
		ReadIncOperation rr_wriop = null;	/* WRITE {@link ReadIncOperation} in rr_interval*/
		
		/** rr_interval: (master_pre_rriop, master_cur_rriop) **/
		for (int rr_index = master_pre_rindex + 1; rr_index < master_cur_rindex; rr_index++)	
		{
			rr_wriop = (ReadIncOperation) master_proc.getOperation(rr_index);
			assertTrue("WRITE ReadIncOperation in rr_interval", rr_wriop.isWriteOp());
			
			rr_wriop.getEarliestRead().initEarlistRead(master_cur_rriop);	// initialize earliest READ
            rr_wriop.setIsPropagated();
			rr_wriop.getLatestWriteMap().updateLatestWrite(pre_riop);	// update latest WRITE map depending on previous operation
			this.riob.getGlobalActiveWritesMap().replace(rr_wriop);	// deactivate some WRITEs
//			System.out.println("Done " + rr_wriop.toString());
			pre_riop = rr_wriop;	// iterate over the next WRITE  
		}
//		System.out.println("ok...");
		master_cur_rriop.getLatestWriteMap().updateLatestWrite(pre_riop);	// update latest WRITE map for @param master_cur_rriop individually
//		System.out.println("ok1....");
		
		// (2) dealing with ww_interval
		ReadIncOperation dw = (ReadIncOperation) master_cur_rriop.getReadfromWrite();	// dictating WRITE for @param master_cur_rriop
		int pid_master = master_cur_rriop.getPid();	// pid of current READ
		int pid_dw = dw.getPid();	// pid of dictating WRITE
		if (pid_dw == pid_master)	// r and D(r) are in the same process and thus ww_interval is empty
		{
		    //如果dw与r在同一线程上，那么相关信息已经更新过了
			if (dw.getIndex() < master_pre_rriop.getIndex())	// D(r) is in r'-downset
				dominated = true;
		}
		else	// r and D(r) are in different processes
		{

		    //dh:把ww interval的概念由po位于D(r)中的操作扩展到D(r)的CausalPast中
            LinkedList<ReadIncOperation> causalPast = causal_set(dw);
            ReadIncOperation causal_w = dw; //防止空值

            Map<String, ReadIncOperation> last_wriop_map = new HashMap<String, ReadIncOperation>();
//            System.out.println("Causal past of " + dw.toString() +":" + causalPast);

            for(int i = 0; i < causalPast.size(); i++){ //按照拓扑序从causalPast中更新
                causal_w = causalPast.get(i);
                update_ww_operation(causal_w, master_cur_rriop, last_wriop_map);
            }
            this.riob.getGlobalActiveWritesMap().addActiveWriteMap(last_wriop_map);
            master_cur_rriop.getLatestWriteMap().updateLatestWrite(causal_w);

            for(ReadIncOperation riop: causalPast){
                riop.resetInCausalSet();
            }

//            if (dw.getIndex() <= dw_pre_wriop.getIndex())	// r and D(r) are in different processes and D(r) is in r'-downset
            if(master_pre_rriop.getPredecessors().contains(dw))
				dominated = true;
//
//        LinkedList<ReadIncOperation> causalPast = causal_set(dw);
//        Map<String, ReadIncOperation> last_wriop_map = new HashMap<String, ReadIncOperation>();
//        System.out.println("Causal past of " + dw.toString() +":" + causalPast);
//        ReadIncOperation causal_w = causalPast.get(0);
//        ReadIncProcess dw_proc = (ReadIncProcess) this.riob.getProcess(causal_w.getPid());// ReadIncProcess in which the first op of dictating WRITE's causal past resides
//        ReadIncOperation dw_pre_wriop = dw_proc.get_pre_wriop();
//        causal_w.getEarliestRead().initEarlistRead(master_cur_rriop);
//        causal_w.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);
//        if(causal_w.isWriteOp() && dw_pre_wriop.isWriteOp()) {
//            this.riob.getGlobalActiveWritesMap().deactivateFrom(dw_pre_wriop, causal_w.getVariable());
//            last_wriop_map.put(causal_w.getVariable(), causal_w);
//        }
//        dw_pre_wriop = causal_w;
//        for(int i = 1; i < causalPast.size(); i++){
//            causal_w = causalPast.get(i);
//            causal_w.getEarliestRead().initEarlistRead(master_cur_rriop);
//            causal_w.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);
//            if(causal_w.isWriteOp()&& dw_pre_wriop.isWriteOp()) {
//                this.riob.getGlobalActiveWritesMap().deactivateFrom(dw_pre_wriop, causal_w.getVariable());
//                last_wriop_map.put(causal_w.getVariable(), causal_w);
//            }
//            dw_pre_wriop = causal_w;
//        }
//        this.riob.getGlobalActiveWritesMap().addActiveWriteMap(last_wriop_map);
//        master_cur_rriop.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);
//
//        if (dw.getIndex() <= dw_pre_wriop.getIndex())	// r and D(r) are in different processes and D(r) is in r'-downset
//            dominated = true;
//        else
//            dw_proc.advance_pre_wriop(dw_pre_wriop);
//


//
//			ReadIncProcess dw_proc = (ReadIncProcess) this.riob.getProcess(pid_dw);	// ReadIncProcess in which dictating WRITE resides
//			ReadIncOperation dw_pre_wriop = dw_proc.get_pre_wriop();	// previous WRITE {@link ReadIncOperation}; changing in iteration
//			assertTrue("Previous ReadIncOperation in ReadIncProcess not with masterPid is WRITE", dw_pre_wriop.isWriteOp());
//
//			ReadIncOperation ww_wriop = null;	// WRITE {@link ReadIncOperation} in ww_interval
//			int dw_index = dw.getIndex();
//			int dw_pre_wriop_index = dw_pre_wriop.getIndex();
//
//			/** dealing with ww_interval: (dw_pre_wriop_index, dw_index] **/
//			Map<String, ReadIncOperation> last_wriop_map = new HashMap<String, ReadIncOperation>();	/** record the last {@link ReadIncOperation} for each variable **/
//			for (int ww_index = dw_pre_wriop_index + 1; ww_index <= dw_index; ww_index++)
//			{
//				ww_wriop = (ReadIncOperation) dw_proc.getOperation(ww_index);
//				System.out.println("ww_wriop:" + ww_wriop.toString());
////				assertTrue("WRITE ReadIncOperation in ww_interval", ww_wriop.isWriteOp()); //dh:因为ww集合现改为w的causalpast，所以也可能出现read操作，故此assert失效
//
//				ww_wriop.getEarliestRead().initEarlistRead(master_cur_rriop);	/** initialize earliest READ **/
//				ww_wriop.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);	/** update LatestWriteMap depending on previous operation **/
//
//				/**
//				 * @description deactivate some WRITEs
//				 * @modified 	hengxin on 2013-1-6
//				 * @reason   	there are two reasons to deactivate WRITEs:
//				 * 			(1) due to dw_pre_wriop (with corresponding READs)
//				 * 			(2) due to ww_interval itself (mainly focus on ones without corresponding READs)
//				 */
//				this.riob.getGlobalActiveWritesMap().deactivateFrom(dw_pre_wriop, ww_wriop.getVariable());	// deactivate some WRITE due to dw_pre_wriop
////                    因为ww_wriop是写操作，所以dw_pre_wriop中对应ww_wriop变量的其他写操作被覆写
////				this.riob.getGlobalActiveWritesMap().deactivateFrom(dw_pre_wriop_constant);  /** deactivate some WRITE due to dw_pre_wriop **/
////				this.riob.getGlobalActiveWritesMap().addActiveWrite(ww_wriop);	/** add this new active WRITE **/
//				last_wriop_map.put(ww_wriop.getVariable(), ww_wriop);	/** record the last {@link ReadIncOperation} for this variable **/
//
//				dw_pre_wriop = ww_wriop;	// iterate over the next WRITE
//			}
//			this.riob.getGlobalActiveWritesMap().addActiveWriteMap(last_wriop_map);	/** add active WRITE for each variable **/
//			master_cur_rriop.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);
//
//			if (dw_index <= dw_pre_wriop_index)	// r and D(r) are in different processes and D(r) is in r'-downset
//				dominated = true;
//			else
//				dw_proc.advance_pre_wriop(dw_pre_wriop);	// advance the previous WRITE forward to the new one


		}
//		/**
//		 *  (3) dealing with @param master_cur_rriop and "dw" separately and specially:
//		 *  	@param master_cur_rriop reads value from "dw", causing other WRITEs with
//		 *  	the same variable are scheduled before "dw" and LatestWrite are updated 
//		 *  	accordingly 
//		 */
//		ReadIncOperation temp_riop = new ReadIncOperation(new RawOperation(GlobalData.WRITE, GlobalData.DUMMYVAR, -1));	// temp
//		for (ReadIncOperation active_wriop : this.riob.getGlobalActiveWritesMap().getActiveWrites(master_cur_rriop.getVariable()))
//			temp_riop.getLatestWriteMap().updateLatestWrite(active_wriop);
////        System.out.println("Finish ww part");
		return dominated;
	}
	
	/**
	 * is it necessary to reschedule operations?
	 * There are two cases in which no reschedule is necessary at all:
	 * (1) D(R) has not been overwritten, and there exist some WRITE other than D(R) 
	 * (2) D(R) has been overwritten, and there still exist some WRITE
	 * 
	 * @return true, if it is necessary to reschedule operations; false, otherwise.
	 * 
	 * @constraints @param dw_cur_rriop has the "writeto" relation with @param cur_rriop
	 * 
	 * @author hengxin
	 * @date 2013-1-9
	 */
	private boolean reschedule_assertion(ReadIncOperation cur_rriop, ReadIncOperation dw_cur_rriop)
	{
		assertTrue("D(R) writes to R", cur_rriop.isReadOp() && dw_cur_rriop.isWriteOp() 
				&& cur_rriop.getVariable().equals(dw_cur_rriop.getVariable())
				&& cur_rriop.getValue() == dw_cur_rriop.getValue());
		
		Set<ReadIncOperation> gaws = this.riob.getGlobalActiveWritesMap().getActiveWrites(cur_rriop.getVariable());
		if(gaws.contains(dw_cur_rriop) && gaws.size() >= 2){
            System.out.println("case 1 cur_r:" + cur_rriop.toString() + " dw:" + dw_cur_rriop.toString());
        }
        else if(! gaws.contains(dw_cur_rriop) && gaws.size() >= 1){
            System.out.println("case 2 cur_r:" + cur_rriop.toString() + " dw:" + dw_cur_rriop.toString());
            System.out.println("gaws:" + gaws.toString());
        }

		return ( (gaws.contains(dw_cur_rriop) && gaws.size() >= 2) || (! gaws.contains(dw_cur_rriop) && gaws.size() >= 1));
	}
	
	/**
	 * @description reschedule (i.e., enforce) orders among {@link ReadIncOperation}s   
	 * 
	 * @param master_cur_rriop current READ {@link ReadIncOperation} being checked			
	 * @return true, if cycle emerges; false, otherwise.
	 */
	private boolean reschedule(ReadIncOperation master_cur_rriop)
	{
		assertTrue("Reschedule operations due to the current READ ReadIncOperation", master_cur_rriop.isReadOp());
		
		ReadIncOperation dw = (ReadIncOperation) master_cur_rriop.getReadfromWrite();	// dictating WRITE {@link ReadIncOperation}
		Set<ReadIncOperation> candidateSet = this.draw_reschedule_boundary(dw);	// identify the possible rescheduled operations
		
		Queue<ReadIncOperation> zeroQueue = new LinkedList<ReadIncOperation>();	// queue containing operations ready to check
		zeroQueue.offer(dw);	// check dw first
		while (! zeroQueue.isEmpty())
		{
			ReadIncOperation wprime_riop = zeroQueue.poll();
			
			/**
			 * @description rule out the case that wprime_riop is READ
			 * 
			 * @modified hengxin on 2013-1-7
			 * @reason wprime_riop maybe READ (e.g., R and D(R) are in the same process)
			 */
			if (wprime_riop.isWriteOp() && ! wprime_riop.isDone())
			{
				// identify the possible W'WR order
				int oldEarlistRead = wprime_riop.getEarliestRead().updateEarliestRead(wprime_riop.getSuccessors());
				ReadIncOperation wriop = wprime_riop.getEarliestRead().identify_wrPair(oldEarlistRead, wprime_riop, this.riob.getMasterProcess());
				
				// apply W'WR order: wprime_riop => wriop
				if (wriop != null)
					if (wprime_riop.apply_wprimew_order(wriop, this.riob))
						return true;	/** cycle **/
				
				// depending on UNDONE operation
				if (wriop != null && wriop.isCandidate() && ! wriop.isDone())
				{
					
//					if (! wriop.isDone())
//					{
						/**
						 * @modified hengxin on 2013-1-9
						 * @modification the update of predecessors and successors 
						 *   are encapsulated in methods related to different kinds of edges
						 * @reason maintenance of predecessors and successors
						 */
//						wprime_riop.getSuccessors().add(wriop);
//						wriop.getPredecessors().add(wprime_riop);
						wprime_riop.incCount();
//					}
//					else
//					{
//						wprime_riop.getEarliestRead().updateEarliestRead(wriop);
//					}
				}
			}
			
			/**
			 * delete dependency and identify operations ready to check
			 * no matter whether it is READ or WRITE 
			 */
			if (wprime_riop.getCount() == 0)
			{
				wprime_riop.setDone();
				
				ReadIncOperation tmp_riop = null;
				for (BasicOperation bop : wprime_riop.getPredecessors())
				{
					tmp_riop = (ReadIncOperation) bop;
					tmp_riop.decCount();
					if (tmp_riop.getCount() == 0 && ! tmp_riop.isDone())	/** TODO is it necessary to check riop.isDone() **/
						zeroQueue.offer(tmp_riop);
				}
			}
		}
		
		/**
		 * @modified hengxin on 2013-1-5
		 */
		for (ReadIncOperation riop : candidateSet)	/** reset flag isCandidate **/
			riop.resetCandidate();
		
		return false;	/** no cycle **/
	}
	
	/**
	 * applying Rule (c): W'WR order due to the fact that 
	 * @param master_cur_rriop must read value from @param dw 
	 * 
	 * @param master_cur_rriop READ {@link ReadIncOperation}
	 * @param dw READ {@link ReadIncOperation}
	 * 
	 * @return true, if cycle emerges; false, otherwise.
	 * 
	 * @modified hengxin on 2013-1-4
	 * @reason you cannot remove elements [apply_wprimew_order] from ArrayList in for each syntax
	 * 
	 * @modified hengxin on 2013-1-9
	 * @modification remove the third parameter {@link ReadIncObservation}
	 * @reason you can assess {@link #riob} in this class directly
	 */
	private boolean readFromDW(ReadIncOperation master_cur_rriop, ReadIncOperation dw)
	{
		assertTrue("READ reads from some WRITE", master_cur_rriop.isReadOp() && dw.isWriteOp());
		
		String var = master_cur_rriop.getVariable();
		
		for (String wriopStr : this.riob.getGlobalActiveWritesMap().getActiveWritesPool(var))
		{
			ReadIncOperation wriop = (ReadIncOperation) this.riob.getWrite(wriopStr);
			if (! wriop.equals(dw) && wriop.apply_wprimew_order(dw,this.riob))	/** apply Rule (c): W'WR order (i.e., wriop => dw => master_cur_rriop) **/
				return true;	/** cycle **/
		}
		
		return false;  /** no cycle **/
	}
	
	/**
	 * @description  perform BFS on G^T from D(r) to identify the 
	 * possible rescheduled {@link ReadIncOperation}s
	 * 
	 * @modified hengxin on 2013-1-5
	 * @reason dealing with {@link ReadIncOperation#isCandidate}
	 * 
	 * @param {@link ReadIncOperation} from which the boundary is drawn
	 * @return set of {@link ReadIncOperation}s involved in the reschedule
	 */
	private Set<ReadIncOperation> draw_reschedule_boundary(ReadIncOperation wriop)
	{
		assertTrue("perform DFS from D(r)", wriop.isWriteOp());
		
		Set<ReadIncOperation> candidateSet = new HashSet<ReadIncOperation>();
		Queue<ReadIncOperation> pending = new LinkedList<ReadIncOperation>();	// pending queue for BFS framework
		wriop.initCount(0);		// the first operation to consider in reschedule (topological sorting) 
		pending.offer(wriop);	// enqueue the start operation
		while (! pending.isEmpty())
		{
			ReadIncOperation cur_op = pending.poll();	// it is a possible rescheduled operation
			cur_op.setCandidate();	// mark the possible rescheduled operation
			candidateSet.add(cur_op);	/** add it to set **/
			cur_op.resetDone();		// reset to be undone
			
			ReadIncOperation tmp_riop = null;
			for (BasicOperation bop: cur_op.getPredecessors())
			{
				tmp_riop = (ReadIncOperation) bop;
				tmp_riop.incCount();	// one more dependency operation
				if (! tmp_riop.isCandidate())
					pending.add(tmp_riop);
			}
			/**
			 *  @modified by hengxin 2013-01-03
			 */
//			pending.addAll(cur_op.getPredecessors());
		}
		
		return candidateSet;
	}

    /**
     *
     * @param wriop
     * @return list of {@link ReadIncOperation}s sorted by their topo order involved in the causalpast of wriop
     */

	private LinkedList<ReadIncOperation> causal_set(ReadIncOperation wriop){
//	    Set<ReadIncOperation> causalSet = new HashSet<>();
//        Queue<ReadIncOperation> pending = new LinkedList<ReadIncOperation>();
//        LinkedList<ReadIncOperation> tempStack = new LinkedList<ReadIncOperation>();
//        LinkedList<ReadIncOperation> opStack = new LinkedList<ReadIncOperation>(); //以逆拓扑序存储的操作集合
//
////        System.out.println("start at wriop:" + wriop.toString());
//
//        pending.offer(wriop);
////        tempStack.addFirst(wriop);
//        while(!pending.isEmpty()){
//            ReadIncOperation curOp = pending.poll();
////            opStack.addFirst(curOp);
////            System.out.println("add " + curOp.toStringCom() + "into tempStack");
//            tempStack.addFirst(curOp); //入栈
//            if(!curOp.isInCausalSet()){ //第一次加入causal set, 初始化入度
//
//                curOp.setInDegree(curOp.getPredecessors().size());
//                curOp.setInCausalSet();
//                causalSet.add(curOp);
////                System.out.println("add " + curOp.toStringCom() +" into causal set, indegree:" + curOp.getInDegree());
//            }
//            ReadIncOperation tempOp = null;
//            for(BasicOperation bop: curOp.getPredecessors()){ //BFS skeleton
//                tempOp = (ReadIncOperation) bop;
////                System.out.println(tempOp.toStringCom() +"is predecessor of " +curOp.toStringCom() +", add it to tempStack");
//                tempStack.addFirst(tempOp);
//                if(! tempOp.isInCausalSet()){
//                    pending.offer(tempOp);
//                }
//            }
//        }
//
//        System.out.println("now finish intialing of causal set");
////        System.out.println("tempStack:" + tempStack.toString());
//
//        for(ReadIncOperation op: tempStack){
//            op.setInCausalSet();
//        }
//
//        while(!tempStack.isEmpty()){
//            ReadIncOperation curOp = tempStack.removeFirst();
////            System.out.println("considering" + curOp.toStringCom()+ " degree:" + curOp.getInDegree() + " isInSet:" + curOp.isInCausalSet());
//            if(curOp.getInDegree() == 0 && curOp.isInCausalSet()){ //找到入度为0且仍在causalset中的点
//                opStack.add(curOp);
//                curOp.resetInCausalSet();
//                ReadIncOperation tempOp = null;
//                for(BasicOperation bop: curOp.getSuccessors()){
//                    tempOp = (ReadIncOperation) bop;
////                    System.out.println("successor of " + curOp.toStringCom() + ":" + tempOp.toStringCom());
//                    if(tempOp.isInCausalSet()){
////                        System.out.println("dec indegree of " + tempOp.toStringCom() + " from " + tempOp.getInDegree());
//                        tempOp.decIndegree();
////                        System.out.println("dec indegree of " + tempOp.toStringCom() + " to " + tempOp.getInDegree());
//                    }
//                    else{
////                        System.out.println(tempOp.toStringCom() + " is not in causal set!");
//                    }
//                }
//            }
//        }
//
//        for(ReadIncOperation op: opStack){ //还原相关信息
////            op.resetInCausalSet();
//            op.setInCausalSet();
//            op.resetIndegree();
//        }
//
//
//        if(opStack.size() != causalSet.size()){
//            System.out.println("opstack size:" + opStack.size() + " causal set size:" + causalSet.size());
//            System.out.println("opstack:" + opStack.toString());
//            System.out.println("causal set:" + causalSet.toString());
//        }
//        else{
//            System.out.println("equal!!!!-------------");
//        }
//

        //restart.
//		System.out.println("start at wriop:" + wriop.toString());
		Set<ReadIncOperation> causalSet = new HashSet<>();
		Queue<ReadIncOperation> pending = new LinkedList<ReadIncOperation>();
		LinkedList<ReadIncOperation> opLi = new LinkedList<ReadIncOperation>(); //以逆拓扑序存储的操作集合
		LinkedList<ReadIncOperation> stack = new LinkedList<>();

		pending.offer(wriop);
		wriop.setInCausalSet();

		ReadIncOperation curOp = null;
		ReadIncOperation tempOp = null;

//		System.out.println("Step1. initial set");
		//第一步，初始化causalset
		while(!pending.isEmpty()) {
			curOp = pending.poll();
			curOp.setInCausalSet();
			causalSet.add(curOp);
			for(BasicOperation bop: curOp.getPredecessors()){ //BFS skeleton
                tempOp = (ReadIncOperation) bop;
//                System.out.println(tempOp.toStringCom() +"is predecessor of " +curOp.toStringCom() +", add it to tempStack");
               	if(!tempOp.isInCausalSet()){
               		pending.offer(tempOp);
				}
            }
		}

//		System.out.println("Step2. init indegree");
		//第二步，初始化入度
		int count = 0;
		for(ReadIncOperation op: causalSet){
			count = 0;
			for(BasicOperation bop: op.getPredecessors()){
				tempOp = (ReadIncOperation) bop;
				if(tempOp.isInCausalSet()){
					count++;
				}
			}
			op.setInDegree(count);
		}

		//找到入度为0的点，入栈，开始拓扑排序
		for(ReadIncOperation op: causalSet){
			if(op.getInDegree() == 0){
				stack.addFirst(op);
				op.setInDegree(-1);
			}
		}

		while (!stack.isEmpty()){
			curOp = stack.removeFirst();
			opLi.add(curOp);
			for(BasicOperation bop: curOp.getSuccessors()){
				tempOp = (ReadIncOperation) bop;
				tempOp.setInDegree(tempOp.getInDegree()-1);
				if(tempOp.getInDegree() == 0){
					stack.addFirst(tempOp);
					tempOp.setInDegree(-1);
				}
			}
		}

		for(ReadIncOperation op: opLi){ //还原相关信息
//            op.resetInCausalSet();
            op.setInCausalSet();
            op.resetIndegree();
        }
//
//		if(opLi.size() != causalSet.size()){
//            System.out.println("opstack size:" + opLi.size() + " causal set size:" + causalSet.size());
//            System.out.println("opstack:" + opLi.toString());
//            System.out.println("causal set:" + causalSet.toString());
//        }
//        else{
//            System.out.println("equal!!!!-------------");
//        }

        assertTrue("Are they equal?", opLi.size() == causalSet.size());
//
	    return opLi;
    }

    public void update_ww_operation(ReadIncOperation curOp, ReadIncOperation master_cur_rriop, Map<String, ReadIncOperation> last_wriop_map){

//        LinkedList<ReadIncOperation> causalPast = causal_set(dw);
//        Map<String, ReadIncOperation> last_wriop_map = new HashMap<String, ReadIncOperation>();
//        System.out.println("Causal past of " + dw.toString() +":" + causalPast);
//        ReadIncOperation causal_w = causalPast.get(0);
//        ReadIncProcess dw_proc = (ReadIncProcess) this.riob.getProcess(causal_w.getPid());// ReadIncProcess in which the first op of dictating WRITE's causal past resides
//        ReadIncOperation dw_pre_wriop = dw_proc.get_pre_wriop();
//        causal_w.getEarliestRead().initEarlistRead(master_cur_rriop);
//        causal_w.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);
//        if(causal_w.isWriteOp() && dw_pre_wriop.isWriteOp()) {
//            this.riob.getGlobalActiveWritesMap().deactivateFrom(dw_pre_wriop, causal_w.getVariable());
//            last_wriop_map.put(causal_w.getVariable(), causal_w);
//        }
//        dw_pre_wriop = causal_w;
//        for(int i = 1; i < causalPast.size(); i++){
//            causal_w = causalPast.get(i);
//            causal_w.getEarliestRead().initEarlistRead(master_cur_rriop);
//            causal_w.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);
//            if(causal_w.isWriteOp()&& dw_pre_wriop.isWriteOp()) {
//                this.riob.getGlobalActiveWritesMap().deactivateFrom(dw_pre_wriop, causal_w.getVariable());
//                last_wriop_map.put(causal_w.getVariable(), causal_w);
//            }
//            dw_pre_wriop = causal_w;
//        }
//        this.riob.getGlobalActiveWritesMap().addActiveWriteMap(last_wriop_map);
//        master_cur_rriop.getLatestWriteMap().updateLatestWrite(dw_pre_wriop);
//
//        if (dw.getIndex() <= dw_pre_wriop.getIndex())	// r and D(r) are in different processes and D(r) is in r'-downset
//            dominated = true;
//        else
//            dw_proc.advance_pre_wriop(dw_pre_wriop);

        if(!curOp.isPropagated()) { //只处理还没有被处理过的写操作 (r-delta中的写)
            curOp.setIsPropagated();
            ReadIncProcess curProcess = (ReadIncProcess) this.riob.getProcess(curOp.getPid());
            ReadIncOperation preOpSameProcess = (ReadIncOperation) curOp.getReProgramOrder();
            if (preOpSameProcess == null) { //如果是某线程第一个操作，就设置为自身
//            System.out.println("preOp itself:" + curOp.toString());
                preOpSameProcess = curOp;
            }
            ReadIncOperation preWriteSameProcess = curProcess.get_pre_wriop();
            curOp.getEarliestRead().initEarlistRead(master_cur_rriop);
            //step1.从同线程前一操作更新latestWriteMap
            curOp.getLatestWriteMap().updateLatestWrite(preOpSameProcess);

            //step2.从causalpast的操作中更新latestWriteMap
            for (BasicOperation visOp : curOp.getPredecessors()) {
                ReadIncOperation visRincOp = (ReadIncOperation) visOp;
                if (visRincOp.isInCausalSet()) {
                    curOp.getLatestWriteMap().updateLatestWrite(visRincOp);
                }
            }

            //step3.根据当前操作更新globalactivewrtiesmap以及当前线程的pre_wriop
            if (curOp.isWriteOp()) {
                if (preOpSameProcess.isWriteOp()) {
                    this.riob.getGlobalActiveWritesMap().deactivateFrom(preOpSameProcess, curOp.getVariable());
                    last_wriop_map.put(curOp.getVariable(), curOp);
                }
                if (curOp.getIndex() > preWriteSameProcess.getIndex()) {
                    curProcess.advance_pre_wriop(curOp);
                }
            }
        }

    }
	
}
