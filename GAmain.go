package main

import (
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/sgal"
	// "code.google.com/p/gosl/mpi"
	// "code.google.com/p/sgal"
	"fmt"
	"scoops3d/Part"
	// IMPORT!!!
)

func main() {

	test := 3000                          //Timeframe
	Pmut, GAout := 0.01, 10 //mutation percentage (1%=0.01)  //GA Results Output Tf // Tolerance
	
	mpi.Start(true)
	// FOS := 5.452500e+00 
	Igood := make([]float64, 5)
	FOS :=  1000.0

	counter := 0.0
	c360 :=0.0
	cimprov:=0.0
	cerr := 0.0

	ovfcn := func(I []float64, arg interface{}) (ov float64, oor int) {

		var (
			score float64
			err   bool
			)

		x, y, z, R, alpha := I[0], I[1], I[2], I[3], I[4]*(180*7)/22
		
		counter = counter + 1.0

		// if alpha > 360.0 {
		// 	score = 10000.0
		// 	c360 = c360+1.0
		// }else{	
		// score, err = Part.Partmain(x, y, z, R, alpha)
		score, err = Part.Partmain(x, y, z, R, alpha, FOS, Igood)
		// fmt.Println(counter)
		// }

		if score<FOS{
			FOS = score
			for i := 0; i < len(Igood); i++ {
				Igood[i]=I[i]				
			}
			fmt.Println(FOS, Igood)
			cimprov = cimprov + 1.0
		}

		if score == 1000.0 {
			cerr = cerr + 1.0
		}
		if err {
			oor = 1 // => out of range
		}

		return score, oor
	}

	report := func(time int, I []float64, arg interface{}) (stop int) {
		fmt.Println("Simulations: ",counter,"FOS: ", FOS, "Improve: ", cimprov, "Over 360: ", c360, "Error: ", cerr)
		return
	}

	var d sgal.Data

	//AriaTagyo
	err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0}, []float64{100.0, 100.0, 100.0, 100.0, 6.2857}) //Creating the grid for possible answers

	//Cone
	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0}, []float64{4000.0, 4000.0, 1500.0, 2500.0, 6.2857}) //Creating the grid for possible answers

	//DonaldGiam
	// err := d.Init([]float64{0.0, 0.0, 0.0, 0.0, 0.0}, []float64{100.0, 100.0, 100.0, 100.0, 6.2857}) //Creating the grid for possible answers

	//Mount StHelen
	// err := d.Init([]float64{557974.631687, 5111446.1291019, 0.0, 0.0, 0.0}, []float64{568974.631687, 5122446.1291019, 4000.0, 5000.0, 360.0}) //Creating the grid for possible answers

	if err != nil {
		fmt.Printf("Failed: %v", err.Error())
		return
	}

	d.Szpop, d.Dtmig, d.Tf, d.Dtout, d.Seltype = 200, 1000, int(test), GAout, 0
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