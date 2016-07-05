// /////////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////RELIABILITY INDEX//////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////
// package main

// import (
// 	// "github.com/cpmech/gosl/mpi"
// 	// "github.com/cpmech/sgal"
// 	"code.google.com/p/gosl/mpi"
// 	"code.google.com/p/sgal"
// 	"fmt"
// 	"math"
// 	"scoops3d/Part"
// 	"time"
// )

// func main() {

// 	//////////////////////////////////////GA INPUT PARAMETER////////////////////////////////////////////////////
// 	test := 200000                          //Timeframe
// 	Pmut, GAout := 0.01, 10 //mutation percentage (1%=0.01)  //GA Results Output Tf // Tolerance
// 	RITolerance, FOSS := 0.01, 1000.0

// 	//////////////////////////////////////MPI STARTS HERE//////////////////////////////////////////////////
// 	mpi.Start(true)
		
// 	//////////////////////////////////////PRINTING TIME////////////////////////////////////////////////////
// 	t := time.Now()
// 	fmt.Println(t)

// 	//////////////////////////////////////MATERIAL PROPERTIES//////////////////////////////////////////////
// 	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()
	
// 	//////////////////////////////////////GA OPTIMISATION PART//////////////////////////////////////////////
// 	ovfcn := func(I []float64, arg interface{}) (ov float64, oor int) {

// 		var (
// 			score float64
// 			score2 float64
// 			err   bool
// 			)
		
// 		///////////////////////////////////////////////////////////////////////////////////////////////////
// 		/////////////////////// CASE STUDY 1 & 3 Mount St Helen and Arai & Tagyo///////////////////////////
// 		///////////////////////////////////////////////////////////////////////////////////////////////////
// 		x, y, z, R, alpha, th1, RI := I[0], I[1], I[2], I[3], I[4], I[5], I[6]
// 		// x, y, z, R, alpha := 1044.6, 2824.2, 1404.8, 1518.2, 141.51
// 		// th1, RI := I[0], I[1]
		

// 		phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)
// 		su11[0] = sukmean[0] + suksd[0]*RI*math.Sin(th1)
		
// 		// ///////////////////////////////////////////////////////////////////////////////////////////////////
// 		// /////////////////////// CASE STDUY 2 DONALD AND GIAM /////////////////////////////////////////////
// 		// ///////////////////////////////////////////////////////////////////////////////////////////////////
// 		// x, y, z, R, alpha, th1, th2, th3, th4, RI := I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8], I[9]

// 		// phi11[2] = phikmean[2] + phiksd[2]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Cos(th4)
// 		// phi11[1] = phikmean[1] + phiksd[1]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Sin(th4)
// 		// phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)*math.Cos(th2)*math.Sin(th3)
// 		// su11[2] = sukmean[2] + suksd[2]*RI*math.Cos(th1)*math.Sin(th2)
// 		// su11[1] = sukmean[1] + suksd[1]*RI*math.Sin(th1)

// 		FOSS, score2, err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)
// 		// FOSS, _, err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)
		
// 		// fmt.Println(I[0], I[1], I[2])

// 		// if math.Abs(FOSS-1.0) < RITolerance && err == false {
// 		if math.Abs(FOSS-1.0) < RITolerance && err == false && score2 > 200000000.0 && score2 < 300000000.0 {		
// 		// if math.Abs(FOSS-1.0) < RITolerance && err == false && score2 > 1000000.0 && score2 < 40000000.0 {					                                                
// 		// if math.Abs(FOSS-1.0) < RITolerance && err == false && score2 > 50000000.0 && score2 < 100000000.0 {		
// 		// if math.Abs(FOSS-1.0) < RITolerance && err == false && score2 > 40000000.0 && score2 < 100000000.0 {
// 			score = RI
// 		} else {
// 			score = RI + 100.0*math.Abs(FOSS-1.0) + 5.0
// 		}
			
// 		if alpha > 360.0{
// 			err = true
// 		}

// 		// if score < 2.6 {
// 		// 	score = 1 + score/10000.0
// 		// }


// 		if err {
// 			oor = 1 // => out of range
// 		}

// 		return score, oor
// 	}

// 	report := func(time int, I []float64, arg interface{}) (stop int) {
// 		// fmt.Println("Mpi:", mpi.Rank())
// 		return
// 	}

// 	var d sgal.Data

// 	// //AriaTagyo
// 	// err := d.Init([]float64{0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 60.0, 360.0, 6.2857, 5.0}) //Creating the grid for possible answers

// 	//Cone
// 	// err := d.Init([]float64{-6000.0, -6000.0, 0.0, 0.0, 0.0, 0.0, 0.0}, []float64{10000.0, 4000.0, 2000.0, 2000.0, 360.0, 6.3, 5.0}) //Creating the grid for possible answers
// 	// err := d.Init([]float64{0.0, 0.0}, []float64{6.3, 5.0}) //Creating the grid for possible answers

// 	// DonaldGiam
// 	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0, -1.5, -1.5, -1.5, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 100.0, 360.0, 1.5, 1.5, 1.5, 6.2857, 5.0}) //Creating the grid for possible answers

// 	//Mount StHelen
// 	// err := d.Init([]float64{557974.631687, 5111446.1291019, 1000.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 180.0}) //Creating the grid for possible answers
// 	err := d.Init([]float64{557974.631687, 5111446.1291019, 1000.0, 0.0, 0.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 360.0, 6.3, 10.0}) //Creating the grid for possible answers
// 	// err := d.Init([]float64{0.0, 0.2}, []float64{6.3, 5.0}) //Creating the grid for possible answers

// 	//Birgham
// 	// err := d.Init([]float64{0.0, 0.0, 1500.0, 0.0, 0.0, 0.0, 0.0}, []float64{3000.0, 3000.0, 3000.0, 3000.0, 360.0, 6.3, 10.0}) //Creating the grid for possible answers

// 	if err != nil {
// 		fmt.Printf("Failed: %v", err.Error())
// 		return
// 	}

// 	d.Szpop, d.Dtmig, d.Tf, d.Dtout, d.Seltype = 200, 10, int(test), GAout, 0
// 	d.Pmut = Pmut

// 	// d.Ranking = true
// 	// err = d.Run(ovfcn, report, nil)

// 	silent := false
// 	err = d.Run(ovfcn, report, silent, nil)
// 	d.UseTime, d.ShowBest, d.Ranking = true, true, true
	

// 	if err != nil {
// 		fmt.Printf("Failed: %v", err.Error())
// 	}
// 	mpi.Stop(true)
	

// }

// 		////////////////////////////////////////////////////////////////////////////
// 		/////////////FACTOR OF SAFETY///////////////////////////////////////////////
// 		////////////////////////////////////////////////////////////////////////////

// package main

// import (
// 	"github.com/cpmech/gosl/mpi"
// 	"github.com/cpmech/sgal"
// 	// "code.google.com/p/gosl/mpi"
// 	// "code.google.com/p/sgal"
// 	"fmt"
// 	"scoops3d/Part"
// 	// "math"
// )

// func main() {

// 	test := 200000                          //Timeframe
// 	Pmut, GAout := 0.01, 10 //mutation percentage (1%=0.01)  //GA Results Output Tf // Tolerance
	
// 	mpi.Start(true)
	
// 	_, _, _, _, gamma, su11, phi11 := Part.SoilP()

// 	ovfcn := func(I []float64, arg interface{}) (ov float64, oor int) {

// 		var (
// 			score float64
// 			// score2 float64
// 			err   bool
// 			)
		

// 		x, y, z, R, alpha := I[0], I[1], I[2], I[3], I[4]

// 		// score, score2,  err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)
// 		score, _,  err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)

// 		if alpha > 360.0{
// 			err = true
// 		}

// 		if score < 1.6 {
// 			score = 1
// 		}


// 		// if score2 > 40000000 || score2 < 20000000{
// 		// 	score = score + math.Abs(30000000-score2)
// 		// }

// 		// if score2 > 10000000 || score2 < 5000000{
// 		// 	score = score + math.Abs(15000000-score2)
// 		// }

// 		// if x > 1600.0 || y > 1000.0 {
// 		// 	score = 100000.0
// 		// }

// 		if err {
// 			oor = 1 // => out of range
// 		}

// 		return score, oor
// 	}

// 	report := func(time int, I []float64, arg interface{}) (stop int) {
// 		fmt.Println("Mpi:", mpi.Rank())
// 		return
// 	}

// 	var d sgal.Data

// 	// //AriaTagyo
// 	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 100.0, 360.0}) //Creating the grid for possible answers
// 	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 100.0, 6.2857, 6.2857, 5.0}) //Creating the grid for possible answers
// 	// err := d.Init([]float64{14.6814, 19.8716, 45.3898, 32.0674, 150.088, 3.09183, 1.72537}, []float64{34.6814, 39.8716, 65.3898, 52.0674, 190.088, 5.09183, 3.72537}) //Creating the grid for possible answers
// 	// err := d.Init([]float64{0.0, 0.0}, []float64{6.3, 10.0}) //Creating the grid for possible answers

// 	//Cone
// 	err := d.Init([]float64{0.0, 0.0, 1200.0, 0.0, 0.0}, []float64{4000.0, 4000.0, 2200.0, 2500.0, 360.0}) //Creating the grid for possible answers

// 	//DonaldGiam
// 	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0}, []float64{100.0, 100.0, 100.0, 100.0, 6.2857}) //Creating the grid for possible answers
// 	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, []float64{100.0, 100.0, 100.0, 100.0, 360.0, 6.3, 1.0}) //Creating the grid for possible answers

// 	//Mount StHelen
// 	// err := d.Init([]float64{557974.631687, 5111446.1291019, 0.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 360}) //Creating the grid for possible answers
	
// 	//Birgham
// 	// err := d.Init([]float64{0.0, 0.0, 1500.0, 0.0, 0.0}, []float64{3000.0, 3000.0, 3000.0, 3000.0, 360.0}) //Creating the grid for possible answers

// 	if err != nil {
// 		fmt.Printf("Failed: %v", err.Error())
// 		return
// 	}

// 	d.Szpop, d.Dtmig, d.Tf, d.Dtout, d.Seltype = 200, 100, int(test), GAout, 0
// 	d.Pmut = Pmut

// 	d.Ranking = true
// 	err = d.Run(ovfcn, report, nil)

// 	// silent := false
// 	// err = d.Run(ovfcn, report, silent, nil)
// 	// d.UseTime, d.ShowBest, d.Ranking = true, true, true
	

// 	if err != nil {
// 		fmt.Printf("Failed: %v", err.Error())
// 	}
// 	mpi.Stop(true)

// }

// 	// ////////////////////////////////////////////////////////////////////////////
// 	// ///////////////////////////////AriaTagyo////////////////////////////////////
// 	// ////////////////////////////////////////////////////////////////////////////
// 	// 1 Material 
// 	// x, y, z, R, alpha :=  24.6543, 34.1117, 50.5277, 37.5, 178.229

// 	// sukmean := []float64{41.65, 41.65, 41.65}
// 	// phikmean:= []float64{15.0, 15.0, 15.0}
// 	// gamma 	:= []float64{18.82, 18.82, 18.82}

// 	// suksd, phiksd := 8.0, 3.0

// 	// su11 := []float64{41.65, 41.65, 41.65}
// 	// phi11:= []float64{15.0, 15.0, 15.0}


// 	// ////////////////////////////////////////////////////////////////////////////
// 	// ///////////////////////////////St.Helens////////////////////////////////////
// 	// ////////////////////////////////////////////////////////////////////////////
// 	// x, y, z, R, alpha :=  561186.0, 5117810.0, 3590.33, 2001.1, 139.631

// /////////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////MARCELO///////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////
// package main

// import (
// 	"github.com/cpmech/gosl/mpi"
// 	"github.com/cpmech/sgal"
// 	// "code.google.com/p/gosl/mpi"
// 	// "code.google.com/p/sgal"
// 	"fmt"
// 	"math"
// 	"scoops3d/Part"
// 	"time"
// )

// func main() {

// 	//////////////////////////////////////GA INPUT PARAMETER////////////////////////////////////////////////////
// 	test := 200000                          //Timeframe
// 	Pmut, GAout := 0.01, 10 //mutation percentage (1%=0.01)  //GA Results Output Tf // Tolerance
// 	RITolerance, FOSS := 0.01, 1000.0

// 	//////////////////////////////////////MPI STARTS HERE//////////////////////////////////////////////////
// 	mpi.Start(true)
		
// 	//////////////////////////////////////PRINTING TIME////////////////////////////////////////////////////
// 	t := time.Now()
// 	fmt.Println(t)

// 	//////////////////////////////////////MATERIAL PROPERTIES//////////////////////////////////////////////
// 	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()
	
// 	//////////////////////////////////////GA OPTIMISATION PART//////////////////////////////////////////////
// 	ovfcn := func(I []float64, arg interface{}) (ov float64, oor int) {

// 		var (
// 			score float64
// 			err   bool
// 			)
		
// 		/////////////////////////////////////////////////////////////////////////////////////////////////////
// 		///////////////////////// CASE STUDY 1 & 3 Mount St Helen and Arai & Tagyo///////////////////////////
// 		/////////////////////////////////////////////////////////////////////////////////////////////////////
// 		th1, RI := I[0], I[1]

// 		x, y, z, R, alpha := 1500.0, 970.0, 2500.0, 550.0, 109.17

// 		phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)
// 		su11[0] = sukmean[0] + suksd[0]*RI*math.Sin(th1)
		
// 		// ///////////////////////////////////////////////////////////////////////////////////////////////////
// 		// /////////////////////// CASE STDUY 2 DONALD AND GIAM /////////////////////////////////////////////
// 		// ///////////////////////////////////////////////////////////////////////////////////////////////////
// 		// x, y, z, R, alpha, th1, th2, th3, th4, RI := I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8], I[9]

// 		// phi11[2] = phikmean[2] + phiksd[2]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Cos(th4)
// 		// phi11[1] = phikmean[1] + phiksd[1]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Sin(th4)
// 		// phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)*math.Cos(th2)*math.Sin(th3)
// 		// su11[2] = sukmean[2] + suksd[2]*RI*math.Cos(th1)*math.Sin(th2)
// 		// su11[1] = sukmean[1] + suksd[1]*RI*math.Sin(th1)

// 		// FOSS, score2, err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)
// 		FOSS, _, err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)
		
// 		if math.Abs(FOSS-1.0) < RITolerance && err == false {
// 		// if math.Abs(FOSS-1.0) < RITolerance && err == false && score2 > 300000000.0 && score2 < 500000000.0 {
// 			score = RI
// 		} else {
// 			// score = RI + 100.0*math.Abs(FOSS-1.0) + 5.0 + 10*math.Abs(score2-10000)
// 			score = RI + 100.0*math.Abs(FOSS-1.0) + 5.0
// 			// fmt.Println(RI, FOSS, score2)
// 		}

// 		// if score2 > 50000 || score2 < 5000{
// 		// 	score = score + math.Abs(25000-score2)
// 		// }
			
// 		if alpha > 360.0{
// 			err = true
// 		}

// 		if err {
// 			oor = 1 // => out of range
// 		}

// 		return score, oor
// 	}

// 	report := func(time int, I []float64, arg interface{}) (stop int) {
// 		fmt.Println("Mpi:", mpi.Rank())
// 		return
// 	}

// 	var d sgal.Data

// 	//AriaTagyo
// 	// err := d.Init([]float64{0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 60.0, 360.0, 6.2857, 5.0}) //Creating the grid for possible answers

// 	//Cone
// 	// err := d.Init([]float64{0.0, 0.0, 1200.0, 0.0, 0.0}, []float64{4000.0, 4000.0, 2200.0, 2500.0, 6.2857}) //Creating the grid for possible answers

// 	//DonaldGiam
// 	// err := d.Init([]float64{0.0, 0.0, 25.0, 0.0, 0.0, -1.5, -1.5, -1.5, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 100.0, 360.0, 1.5, 1.5, 1.5, 6.2857, 5.0}) //Creating the grid for possible answers

// 	//Mount StHelen
// 	// err := d.Init([]float64{557974.631687, 5111446.1291019, 1000.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 180.0}) //Creating the grid for possible answers
// 	// err := d.Init([]float64{557974.631687, 5111446.1291019, 1000.0, 0.0, 0.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 360.0, 6.3, 10.0}) //Creating the grid for possible answers
// 	err := d.Init([]float64{0.0, 0.0}, []float64{6.3, 5.0}) //Creating the grid for possible answers

// 	if err != nil {
// 		fmt.Printf("Failed: %v", err.Error())
// 		return
// 	}

// 	d.Szpop, d.Dtmig, d.Tf, d.Dtout, d.Seltype = 200, 10, int(test), GAout, 0
// 	d.Pmut = Pmut

// 	d.Ranking = true
// 	err = d.Run(ovfcn, report, nil)

// 	// silent := false
// 	// err = d.Run(ovfcn, report, silent, nil)
// 	// d.UseTime, d.ShowBest, d.Ranking = true, true, true
	

// 	if err != nil {
// 		fmt.Printf("Failed: %v", err.Error())
// 	}
// 	mpi.Stop(true)
	

// }

/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////FOS//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
package main

import (
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/sgal"
	// "code.google.com/p/gosl/mpi"
	// "code.google.com/p/sgal"
	"fmt"
	// "math"
	"scoops3d/Part"
	"time"
)

func main() {

	//////////////////////////////////////GA INPUT PARAMETER////////////////////////////////////////////////////
	test := 200000                          //Timeframe
	Pmut, GAout := 0.01, 10 //mutation percentage (1%=0.01)  //GA Results Output Tf // Tolerance
	FOSS := 1000.0

	//////////////////////////////////////MPI STARTS HERE//////////////////////////////////////////////////
	mpi.Start(true)
		
	//////////////////////////////////////PRINTING TIME////////////////////////////////////////////////////
	t := time.Now()
	fmt.Println(t)

	//////////////////////////////////////MATERIAL PROPERTIES//////////////////////////////////////////////
	_, _, _, _, gamma, su11, phi11 := Part.SoilP()
	
	//////////////////////////////////////GA OPTIMISATION PART//////////////////////////////////////////////
	ovfcn := func(I []float64, arg interface{}) (ov float64, oor int) {

		var (
			score float64
			// score2 float64
			err   bool
			)
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////// CASE STUDY 1 & 3 Mount St Helen and Arai & Tagyo///////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////
		x, y, z, R, alpha, nc1, nc2, nc3, nc4, nc5 := I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8], I[9]		

		// nc1, nc2, nc3, nc4, nc5	= 0.0, 0.0, 0.0, 0.0, 0.0
		// ///////////////////////////////////////////////////////////////////////////////////////////////////
		// /////////////////////// CASE STDUY 2 DONALD AND GIAM /////////////////////////////////////////////
		// ///////////////////////////////////////////////////////////////////////////////////////////////////
		// x, y, z, R, alpha, th1, th2, th3, th4, RI := I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8], I[9]

		// phi11[2] = phikmean[2] + phiksd[2]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Cos(th4)
		// phi11[1] = phikmean[1] + phiksd[1]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Sin(th4)
		// phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)*math.Cos(th2)*math.Sin(th3)
		// su11[2] = sukmean[2] + suksd[2]*RI*math.Cos(th1)*math.Sin(th2)
		// su11[1] = sukmean[1] + suksd[1]*RI*math.Sin(th1)

		FOSS, _, err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma, nc1, nc2, nc3, nc4, nc5)
		// FOSS, _, err = Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma, nc1, nc2, nc3, nc4, nc5)
		
		if err == false {		
			score = FOSS
		} else {
			score = 10000.0
		}
			
		if alpha > 360.0{
			err = true
		}

		if err {
			oor = 1 // => out of range
		}

		return score, oor
	}

	report := func(time int, I []float64, arg interface{}) (stop int) {
		// fmt.Println("Mpi:", mpi.Rank())
		return
	}

	var d sgal.Data

	// //AriaTagyo
	err := d.Init([]float64{0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 60.0, 360.0, 5.0, 5.0, 5.0, 5.0, 5.0}) //Creating the grid for possible answers

	//Cone
	// err := d.Init([]float64{-6000.0, -6000.0, 0.0, 0.0, 0.0, 0.0, 0.0}, []float64{10000.0, 4000.0, 2000.0, 2000.0, 360.0, 6.3, 5.0}) //Creating the grid for possible answers
	// err := d.Init([]float64{0.0, 0.0}, []float64{6.3, 5.0}) //Creating the grid for possible answers

	// DonaldGiam
	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0, -1.5, -1.5, -1.5, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 100.0, 360.0, 1.5, 1.5, 1.5, 6.2857, 5.0}) //Creating the grid for possible answers

	//Mount StHelen
	// err := d.Init([]float64{557974.631687, 5111446.1291019, 1000.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 180.0}) //Creating the grid for possible answers
	// err := d.Init([]float64{557974.631687, 5111446.1291019, 1000.0, 0.0, 0.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 360.0, 6.3, 10.0}) //Creating the grid for possible answers
	// err := d.Init([]float64{0.0, 0.2}, []float64{6.3, 5.0}) //Creating the grid for possible answers

	//Birgham
	// err := d.Init([]float64{0.0, 0.0, 1500.0, 0.0, 0.0, 0.0, 0.0}, []float64{3000.0, 3000.0, 3000.0, 3000.0, 360.0, 6.3, 10.0}) //Creating the grid for possible answers

	if err != nil {
		fmt.Printf("Failed: %v", err.Error())
		return
	}

	d.Szpop, d.Dtmig, d.Tf, d.Dtout, d.Seltype = 200, 10, int(test), GAout, 0
	d.Pmut = Pmut

	d.Ranking = true
	err = d.Run(ovfcn, report, nil)

	// silent := false
	// err = d.Run(ovfcn, report, silent, nil)
	// d.UseTime, d.ShowBest, d.Ranking = true, true, true
	

	if err != nil {
		fmt.Printf("Failed: %v", err.Error())
	}
	mpi.Stop(true)
}
