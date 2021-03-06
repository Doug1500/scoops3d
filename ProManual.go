// /////////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////RELIABILITY INDEX//////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////
// package main

// import (
// 	"fmt"
// 	"scoops3d/Part"
// 	"math"
// )

// func main() {

// 	fmt.Println("Welcome to GASlope 3D LEM")
// 	fmt.Println("Can you please x, y, z, R, alpha, th1, RI")

// 	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()

// 	RITolerance := 0.01

// 	// /////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// ///////////////////////// CASE STUDY 1 & 3 Mount St Helen and Arai & Tagyo///////////////////////////
// 	// /////////////////////////////////////////////////////////////////////////////////////////////////////
// 	// var x, y, z, R, alpha, th1, RI float64
// 	// fmt.Scan(&x, &y, &z, &R, &alpha, &th1, &RI)

// 	// // x, y, z, R, alpha := 1500.0, 970.0, 2500.0, 550.0, 109.17
// 	// // var th1, RI float64
// 	// // fmt.Scan(&th1, &RI)
	
// 	// phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)
// 	// su11[0] = sukmean[0] + suksd[0]*RI*math.Sin(th1)

// 	///////////////////////////////////////////////////////////////////////////////////////////////////
// 	/////////////////////// CASE STDUY 2 DONALD AND GIAM /////////////////////////////////////////////
// 	///////////////////////////////////////////////////////////////////////////////////////////////////
// 	var x, y, z, R, alpha, th1, th2, th3, th4, RI float64
// 	fmt.Scan(&x, &y, &z, &R, &alpha, &th1, &th2, &th3, &th4, &RI)
	
// 	phi11[2] = phikmean[2] + phiksd[2]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Cos(th4)
// 	phi11[1] = phikmean[1] + phiksd[1]*RI*math.Cos(th1)*math.Cos(th2)*math.Cos(th3)*math.Sin(th4)
// 	phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)*math.Cos(th2)*math.Sin(th3)
// 	su11[2] = sukmean[2] + suksd[2]*RI*math.Cos(th1)*math.Sin(th2)
// 	su11[1] = sukmean[1] + suksd[1]*RI*math.Sin(th1)

		
// 	FOSS, _, err := Part.Partmain(x, y, z, R, alpha, 0, su11, phi11, gamma)

// 	score := 10000.0
// 	if math.Abs(FOSS-1.0) < RITolerance && err == false {
// 		score = RI
// 	} else {
// 		score = RI + 100.0*math.Abs(FOSS-1.0) + 5.0
// 	}

// 	fmt.Println(phi11, su11, sukmean)
// 	fmt.Println(FOSS)
// 	fmt.Println(score, err)
// 	fmt.Println("RI:", RI)

// }


// // /////////////////////////////////////////////////////////////////////////////////////
// // //////////////////////FACTOR OF SAFETY///////////////////////////////////////////////
// // ////////////////////////////////////////////////////////////////////////////////////
// // package main

// // import (
// // 	"fmt"
// // 	"scoops3d/Part"
// // )

// // func main() {

// // 	fmt.Println("Welcome to GASlope 3D LEM")
// // 	fmt.Println("Can you please x, y, z, R, alpha, th1, RI")

// // 	_, _, _, _, gamma, su11, phi11 := Part.SoilP()

// // 	var x, y, z, R, alpha float64
// // 	fmt.Scan(&x, &y, &z, &R, &alpha)

// // 	FOSS,_, err := Part.Partmain(x, y, z, R, alpha, 0, su11, phi11, gamma)
	
	
// // 	fmt.Println(FOSS, err)
	

// // }


/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////FOS//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
package main

import (
	"fmt"
	"math"
	"scoops3d/Part"
	"time"
)

func main() {

	//////////////////////////////////////PRINTING TIME////////////////////////////////////////////////////
	t := time.Now()
	fmt.Println(t)

	//////////////////////////////////////MATERIAL PROPERTIES//////////////////////////////////////////////
	_, _, _, _, gamma, su11, phi11 := Part.SoilP()

	var x, y, z, R, alpha, nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10 float64
	fmt.Scan(&x, &y, &z, &R, &alpha, &nc1, &nc2, &nc3, &nc4, &nc5, &nc6, &nc7, &nc8, &nc9, &nc10)	
	nc1 = 0.0
	// x, y, z, R, alpha := 26.9625, 49.65, 54.395, 32.7162, 180.0
	// nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10	:= 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	// nc1, nc2, nc3, nc4, nc5	:= 0.0, 0.0, 7.0, 5.0, 20.0
	// nc1, nc2, nc3, nc4, nc5	:= 0.0, 0.566502, 1.08243, 1.82565, 2.55044
	// nc1, nc2, nc3, nc4, nc5	:= 10.0, 0.66712, 8.01744, 0.4856, 9.70392
	// x, y, z, R, alpha = 26.9625, 49.65, 54.395, 32.7162, 181.084
	I := []float64{x, y, z, R, alpha, nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10}
	angflag := false
		for i := 5; i < len(I)-1; i++ {
			fmt.Println(math.Abs(I[i+1]-I[i]))
			if math.Abs(I[i+1]-I[i]) > 5.0 || math.Abs(I[i+1]-I[i]) < 0.0 {
				angflag = true
				i = len(I)
			}			
		}

	// FOSS, score2, err = Part.Partmain(x, y, z, R, alpha, mpi.Rank(), su11, phi11, gamma)
	score, _, err := Part.Partmain(x, y, z, R, alpha, 16, su11, phi11, gamma, nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10)
	
	fmt.Println(score, err,angflag)

}
