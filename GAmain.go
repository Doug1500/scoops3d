package main

import (
	// "github.com/cpmech/gosl/mpi"
	// "github.com/cpmech/sgal"
	"code.google.com/p/gosl/mpi"
	"code.google.com/p/sgal"
	"fmt"
	"math"
	"scoops3d/Part"
)

func main() {

	test := 200000                          //Timeframe
	Pmut, GAout := 0.01, 10 //mutation percentage (1%=0.01)  //GA Results Output Tf // Tolerance
	
	mpi.Start(true)

	RITolerance := 0.01
	FOSS := 1000.0

	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()
	// _, _, _, _, gamma, su11, phi11 := Part.SoilP()

	ovfcn := func(I []float64, arg interface{}) (ov float64, oor int) {

		var (
			score float64
			err   bool
			)
		
		//////////////////////////////////////////////////////////////////////////////
		///////////////FACTOR OF SAFETY///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////
		// x, y, z, R, alpha := I[0], I[1], I[2], I[3], I[4]*(180*7)/22

		//////////////////////////////////////////////////////////////////////////////
		///////////////RELIABILITY INDEX//////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////
		x, y, z, R, alpha, th1, RI := I[0], I[1], I[2], I[3], I[4]*180*7/22, I[5], I[6]

		phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)
		su11[0] = sukmean[0] + suksd[0]*RI*math.Sin(th1)
	
		FOSS, err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)
		// fmt.Println(phi11[1], su11[1], FOSS)
		if math.Abs(FOSS-1.0) < RITolerance && err == false {
			score = RI
		} else {
			score = RI + 100.0*math.Abs(FOSS-1.0) + 5.0
		}


		//////////////////////////////////////////////////////////////////////////////
		///////////////FACTOR OF SAFETY///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////
		// score, err = Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma)

		if err {
			oor = 1 // => out of range
		}

		return score, oor
	}

	report := func(time int, I []float64, arg interface{}) (stop int) {
		fmt.Println("Mpi:", mpi.Rank())
		return
	}

	var d sgal.Data

	//AriaTagyo
	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 100.0, 6.2857}) //Creating the grid for possible answers
	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, []float64{60.0, 60.0, 100.0, 100.0, 6.2857, 6.2857, 5.0}) //Creating the grid for possible answers
	// err := d.Init([]float64{14.6814, 19.8716, 45.3898, 32.0674, 150.088, 3.09183, 1.72537}, []float64{34.6814, 39.8716, 65.3898, 52.0674, 190.088, 5.09183, 3.72537}) //Creating the grid for possible answers
	// err := d.Init([]float64{0.0, 0.0}, []float64{6.3, 10.0}) //Creating the grid for possible answers

	//Cone
	// err := d.Init([]float64{0.0, 0.0, 1200.0, 0.0, 0.0}, []float64{4000.0, 4000.0, 2200.0, 2500.0, 6.2857}) //Creating the grid for possible answers

	//DonaldGiam
	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0}, []float64{100.0, 100.0, 100.0, 100.0, 6.2857}) //Creating the grid for possible answers
	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, []float64{100.0, 100.0, 100.0, 100.0, 360.0, 6.3, 1.0}) //Creating the grid for possible answers

	//Mount StHelen
	// err := d.Init([]float64{557974.631687, 5111446.1291019, 0.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 180.0}) //Creating the grid for possible answers
	err := d.Init([]float64{557974.631687, 5111446.1291019, 0.0, 0.0, 0.0, 0.0, 0.2}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 180.0, 6.3, 10.0}) //Creating the grid for possible answers
	// err := d.Init([]float64{0.0, 0.2}, []float64{6.3, 5.0}) //Creating the grid for possible answers

	if err != nil {
		fmt.Printf("Failed: %v", err.Error())
		return
	}

	d.Szpop, d.Dtmig, d.Tf, d.Dtout, d.Seltype = 200, 10, int(test), GAout, 0
	d.Pmut = Pmut

	// d.Ranking = true
	// err = d.Run(ovfcn, report, nil)

	silent := false
	err = d.Run(ovfcn, report, silent, nil)
	d.UseTime, d.ShowBest, d.Ranking = true, true, true
	

	if err != nil {
		fmt.Printf("Failed: %v", err.Error())
	}
	mpi.Stop(true)

}
	//////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////AriaTagyo////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	//1 Material 
	// x, y, z, R, alpha :=  24.6543, 34.1117, 50.5277, 37.5, 178.229

	// sukmean := []float64{41.65, 41.65, 41.65}
	// phikmean:= []float64{15.0, 15.0, 15.0}
	// gamma 	:= []float64{18.82, 18.82, 18.82}

	// suksd, phiksd := 8.0, 3.0

	// su11 := []float64{41.65, 41.65, 41.65}
	// phi11:= []float64{15.0, 15.0, 15.0}


	//////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////St.Helens////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	// x, y, z, R, alpha :=  561186.0, 5117810.0, 3590.33, 2001.1, 139.631
