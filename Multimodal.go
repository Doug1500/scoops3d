package main

import (
	// "github.com/cpmech/gosl/mpi"
	// "github.com/cpmech/sgal"
	"code.google.com/p/gosl/mpi"
	"code.google.com/p/sgal"
	"fmt"
	"math"
)

func main() {

	// Trials := 10
	// MPRes := make([]float64, Trials)
	MPRes := 0.0

	mpi.Start(true)

	// for i := 0; i < Trials; i++ {
	ovfcn := func(I []float64, arg interface{}) (ov float64, oor int) {

		var (
			score float64
			err bool
			)

		x := I[0]

		// if x >= 1.0 || x<=0.0{
		// 	err = true
		// }

		// C := math.Sin(2.0*(22.0/7.0)*x)*math.Sin(2.0*(22.0/7.0)*x)
		// B := math.Exp(-x*x)
		// A := (B)*(C)
		// score = 1.0 / (1.0 - A)



		C := math.Sin(5.0*(22.0/7.0)*x)*math.Sin(5.0*(22.0/7.0)*x)
		B := math.Exp(-x*x)
		A := (B)*(C)
		score = 1.0 / (1.1 - A)

		if x >= 0.9 || x <= 0.05{
			// err = true
			score = 1000.0
		}

		if err {
			oor = 1 // => out of range
		}

	return score, oor
	}

	report := func(time int, I []float64, arg interface{}) (stop int) {
		if time == 200 {
			// MPRes[i] = I[0]
			MPRes = I[0]
		}
		return
	}

	var d sgal.Data

	err := d.Init([]float64{-10}, []float64{10}) //Creating the grid for possible answers

	if err != nil {
		fmt.Printf("Failed: %v", err.Error())
		return
	}

	d.Szpop, d.Dtmig, d.Tf, d.Dtout, d.Seltype = 10, 10, 1000, 10, 0
	d.Pmut = 0.01

	// d.Ranking = true
	// err = d.Run(ovfcn, report, nil)

	silent := false
	err = d.Run(ovfcn, report, silent, nil)
	d.UseTime, d.ShowBest, d.Ranking = true, true, true
	
	if err != nil {
		fmt.Printf("Failed: %v", err.Error())
	}
			
	// }
	mpi.Stop(true)
	fmt.Println(MPRes)
	

}