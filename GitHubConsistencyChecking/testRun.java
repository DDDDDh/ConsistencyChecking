import src.cn.edu.nju.moon.consistency.model.*;
import src.cn.edu.nju.moon.consistency.checker.*;
import src.cn.edu.nju.moon.consistency.datastructure.*;
import src.cn.edu.nju.moon.consistency.schedule.*;
import src.cn.edu.nju.moon.consistency.model.observation.constructor.*;
import src.cn.edu.nju.moon.consistency.model.observation.*;
import src.cn.edu.nju.moon.consistency.schedule.constructor.*;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import java.util.*;



public class testRun {


    public static void main (String args[]){

        System.out.println("Begin");

        GlobalData.VISUALIZATION = false;

//        IBasicObservationConstructor randcons = new FileBasicObservationConstructor("/Users/yi-huang/Project/IncrementalCCC/target/ParameterChoosing/Running_202151410_opNum1000_processNum10_varRange20_valRange100_rRate3_wRate1_debug.edn");
        IBasicObservationConstructor randcons = new FileBasicObservationConstructor("/Users/yi-huang/Project/IncrementalCCC/target/ParameterChoosing/debug/special01_debug.edn");
        BasicObservation bob = randcons.construct();


        long startTime = System.nanoTime();

        Checker cl_checker_rand = new ClosureGraphChecker(bob, randcons.get_ob_id() + "check", new WeakSchedule(bob.getProcNum()));
        System.out.println("here");
        boolean isPRAM = cl_checker_rand.check();
        System.out.println("Satisfied PRAM:" + isPRAM);

        long endTime = System.nanoTime();

        System.out.println("Checking time for b-f alg: " + (endTime-startTime) );

//        assertFalse("Random 1459 does not satisfy PRAM Consistency", cl_checker_rand.check());
//        assertTrue("Process 0 and Process 1 has legal view for PRAM Consistency",
//                Arrays.equals(((WeakSchedule) cl_checker_rand.getSchedule()).constructBinarySchedule(),
//                        new Boolean[] {true, false, false, false, false, false, false, false, false, false}));
//        assertFalse("The views in schedule are all valid", cl_checker_rand.getSchedule().valid());
        System.out.println(cl_checker_rand.getSchedule());

//        System.out.println("Test generate random view");
//        RandomValidViewFactory randomValidViewFactory = new RandomValidViewFactory();
//        View tempView = randomValidViewFactory.generateView(5,10,50);
//        System.out.println(tempView.toString());

        startTime = System.nanoTime();
        Checker incre_checker = new ReadIncChecker(bob, randcons.get_ob_id() + "check", new WeakSchedule(bob.getProcNum()));
        isPRAM = incre_checker.check();
        System.out.println("Satisfied PRAM 2:" + isPRAM);
        endTime = System.nanoTime();
        System.out.println("Checking time for inc alg: " + (endTime-startTime));

    }




}
