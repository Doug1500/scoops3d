package Part

import (
	// "fmt"
	"math"
	"io/ioutil"
	"strconv"
	"os/exec"
	"bufio"
	"log"
	"strings"
	"os"
)

func Partmain(x, y, z, R, alpha float64, mprank int, su, phi, gamma []float64, nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10 float64) (score, score2 float64, err bool){
// func Partmain(x, y, z, R, alpha float64,rank int, su, phi, gamma []float64 ) (score float64, err bool){
// func Partmain(x, y, z, R, alpha, FOS float64, Igood []float64) (score float64, err bool){
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////SINGLE SOIL MATERIAL/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	rank :=	strconv.Itoa(mprank)
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material := "1   41.65   15   18.82"
			 // tag	Cohesion phi gamma	
	

	// FileName := "StHelens"
	// ASCFile := "sthel_res100mDEM.asc"
	// material := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material := "1   1000.0   40.0   24.0"
	
	// FileName := "Cone"
	// ASCFile := "cone1000DEM.asc"
	// material := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material := "1   17630   40   21"
	// material := "1   1444   30   20"
	// material := "1   172.9   20   19"

	// FileName := "DonaldGiam"
	// ASCFile := "emb10DEM.asc"
	// material := "1   3   19.6   20"

	// FileName := "Birgham"
	// ASCFile := "emb20DEM.asc"
	// material := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +" " + strconv.FormatFloat(nc1,'E',-1, 64) +" "+ strconv.FormatFloat(nc2,'E',-1, 64)+" "+ strconv.FormatFloat(nc3,'E',-1, 64) +" "+ strconv.FormatFloat(nc4,'E',-1, 64) +" "+ strconv.FormatFloat(nc5,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\noutput directory\nOutput/")

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////MULTIPLE MATERIAL/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	FileName := "AraiTagyo"
	ASCFile := "emb20DEM.asc"
	LayerFile := "emb20layer"	
	material1 := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	material2 := "2   "+strconv.FormatFloat(su[1],'E',-1, 64)+"   "+strconv.FormatFloat(phi[1],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[1],'E',-1, 64)
	material3 := "3   "+strconv.FormatFloat(su[2],'E',-1, 64)+"   "+strconv.FormatFloat(phi[2],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[2],'E',-1, 64)

	// material1 := "1   29.4   12   18.82"
	// material2 := "2   9.8   5   18.82"
	// material3 := "3   294   40   18.82" 


	// FileName := "DonaldGiam"
	// ASCFile := "emb10DEM.asc"
	// LayerFile := "emb10layer"
	// material1 := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material2 := "2   "+strconv.FormatFloat(su[1],'E',-1, 64)+"   "+strconv.FormatFloat(phi[1],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[1],'E',-1, 64)
	// material3 := "3   "+strconv.FormatFloat(su[2],'E',-1, 64)+"   "+strconv.FormatFloat(phi[2],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[2],'E',-1, 64)

	// material1 := "1   0   38   19.5"
	// material2 := "2   5.3   23   19.5"
	// material3 := "3   7.2   20   19.5" 
	// LayerFile := "emb10layer"

	d1 := []byte("title\nScoops3D example R; Mount Saint StHelensens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nnmat\n3\nlnum   cee   phi   gamt\n"+material1+"\n"+material2+"\n"+material3+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +" "+ strconv.FormatFloat(nc1,'E',-1, 64) +" "+ strconv.FormatFloat(nc2,'E',-1, 64)+" "+ strconv.FormatFloat(nc3,'E',-1, 64) +" "+ strconv.FormatFloat(nc4,'E',-1, 64) +" "+ strconv.FormatFloat(nc5,'E',-1, 64) +" "+ strconv.FormatFloat(nc6,'E',-1, 64) +" "+ strconv.FormatFloat(nc7,'E',-1, 64) +" "+ strconv.FormatFloat(nc8,'E',-1, 64) +" "+ strconv.FormatFloat(nc9,'E',-1, 64) +" "+ strconv.FormatFloat(nc10,'E',-1, 64) +" 1000.0\nremove   foscut\nM   2.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\nlayer file\n"+FileName+"/input/"+LayerFile+"\noutput directory\nOutput/")
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////3D MULTIPLE MATERIAL DATA////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   -1   -1   -1"
	// MaterialFile := "emb20mat3D.txt"

	// FileName := "DonaldGiam"
	// ASCFile := "emb10DEM.asc"
	// material := "1   -1   -1   -1"
	// MaterialFile := "emb10mat3d.txt"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nstr3d   linterp\n1   0\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   5.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\nmaterial properties file\n"+FileName+"/input/"+MaterialFile+"\noutput directory\nOutput/")	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////PIEZOMETER HEAD DATA/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   41.65   15   18.82   18.82"
	// PiezoFile := "emb20piezo.asc"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater  gamw\npz   9.81\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\npiezometric file\n"+FileName+"/input/"+PiezoFile+"\noutput directory\nOutput/")	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////PRESSURE HEAD DATA///////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   41.65   15   18.82   18.82"
	// Piezo3DFile := "emb20phead3D.txt"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater  gamw\n3d   9.81\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\npressure head file\n"+FileName+"/input/"+Piezo3DFile+"\noutput directory\nOutput/")	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////MATERIAL DATA RU///////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FileName := "Cone"
	// ASCFile := "cone1000DEM.asc"
	// material := "1   1444   30   20   0.2"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nru\nnmat\n1\nlnum   cee   phi   gamt   ru\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n1\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\noutput directory\nOutput/")
	err1 := ioutil.WriteFile("test"+rank+".scp", d1, 0644)
	// err1 := ioutil.WriteFile("test.scp", d1, 0644)
	Check(err1)

	//exe.go --> 
	exec.Command("/home/yewintun/mygo/src/scoops3d/fortrancode/scoops3d"+rank+".exe").Run()
	
	//read.go --> read test_out.txt -->FOS
	file, err3 := os.Open("/home/yewintun/mygo/src/scoops3d/Output/test"+rank+"_out.txt")
	// file, err3 := os.Open("/home/yewintun/mygo/src/scoops3d/Output/test_out.txt")

	if err3 != nil {
		log.Fatal(err3)
	}

	defer file.Close()

	important, important2, important4 := "1st", "2nd", "3rd"
	counter := 0.0

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		a:=scanner.Text()
		if counter < 2.0{
			if a == "3D POTENTIAL FAILURE"{
				important = a
			}
			if important == "3D POTENTIAL FAILURE"{
				important2 = a
				counter = counter + 1.0
			}
		}else{
			counter = counter + 1.0
			if counter == 4.0{
				important4 = a
			}
		}
	}
	
	important3 := strings.Fields(important2)
	important5 := strings.Fields(important4)
	
	
	if important3[0] == "2nd" {
		score = 1000.0
	}else{
		score, _ =(strconv.ParseFloat(important3[5], 64))
		if score == 0.0{
			score = 1000.0
		}
	}

	if important5[0] == "3rd" {
		score2 = 1000.0
	}else{
		score2, _ =(strconv.ParseFloat(important5[3], 64))
		// fmt.Println(score2)
		if score2 == 0.0{
			score2 = 10.0
		}
	}


	if err3 := scanner.Err(); err3 != nil {
		log.Fatal(err3)
	}

	if err3 := scanner.Err(); err3 != nil {
		log.Fatal(err3)
	}

	// if score== 1000.0 || score == 100.0{
	// 	I := []float64{x, y, z, R, alpha}
	// 	score = SetValue(FOS, Igood, I)
	// }

	return score, score2, err
}

func Check(e error) {
	if e != nil {
		panic(e)
	}
}

func SetValue(FOS float64, Igood, I []float64) float64 {
	Ifinal := make([]float64, len(Igood))
	for i := 0; i < len(Igood); i++ {
		Ifinal[i] = (Igood[i]-I[i])*(Igood[i]-I[i])
	}
	score := math.Sqrt(SUM(Ifinal))*1000.0 + FOS
	// fmt.Println(score, I)
	return score
}

func SUM(slice []float64) float64 {
	sum := 0.0
	for i := 0; i < len(slice); i++ {
		sum = slice[i] + sum
	}
	return sum
}

func SoilP() ([]float64, []float64, []float64, []float64, []float64, []float64, []float64) {
	// // // ////////////////////////////////////////////////////////////////////////////
	// // // ///////////////////////////////AriaTagyo////////////////////////////////////
	// // // ////////////////////////////////////////////////////////////////////////////
	// sukmean := []float64{41.65, 41.65, 41.65}
	// phikmean:= []float64{15.0, 15.0, 15.0}
	// gamma 	:= []float64{18.82, 18.82, 18.82}

	// suksd := []float64{8.0, 0.0, 0.0}
	// phiksd:= []float64{3.0, 0.0, 0.0}
	
	// su11 := []float64{41.65, 41.65, 41.65}
	// phi11:= []float64{15.0, 15.0, 15.0}

	//////////////////////////////////////////////////////////////////////////
	/////////////////////////////AT Multiplelayer//////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	sukmean := []float64{29.4, 9.8, 294.0}
	phikmean:= []float64{12.0, 5.0, 40.0}
	

	suksd := []float64{29.4, 9.8, 294.0}
	phiksd:= []float64{12.0, 5.0, 40.0}
	
	gamma 	:= []float64{18.82, 18.82, 18.82}
	su11 := []float64{20.8, 400.8, 400.0}
	phi11:= []float64{18.0, 30.0, 30.0}

	// gamma 	:= []float64{18.82, 18.82, 18.82}
	// su11 := []float64{20.8, 20.8, 400.0}
	// phi11:= []float64{18.0, 18.0, 30.0}

	// //////////////////////////////////////////////////////////////////////////
	// /////////////////////////////Donald&Giam//////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////
	// sukmean := []float64{0.0, 5.3, 7.2}
	// phikmean:= []float64{38.0, 23.0, 20.0}
	// gamma 	:= []float64{19.5, 19.5, 19.5}

	// suksd := []float64{0.0, 0.53, 1.44}
	// phiksd:= []float64{3.8, 4.6, 4.0}
	
	// su11 := []float64{0.0, 5.3, 7.2}
	// phi11:= []float64{38.0, 23.0, 20.0}

	////////////////////////////////////////////////////////////////////////////
	///////////////////////////////St.Helens////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	// sukmean := []float64{1000.0, 1000.0, 1000.0}
	// phikmean:= []float64{40.0, 40.0, 40.0}
	// gamma 	:= []float64{24.0, 24.0, 24.0}

	// suksd := []float64{230.0, 0.0, 0.0}
	// phiksd:= []float64{2.8, 0.0, 0.0}	
	
	// su11 := []float64{1000.0, 1000.0, 1000.0}
	// phi11:= []float64{40.0, 40.0, 40.0}

	// //////////////////////////////////////////////////////////////////////////
	// /////////////////////////////Birgham//////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////
	// sukmean := []float64{590.0, 590.0, 590.0}
	// phikmean:= []float64{35.0, 35.0, 35.0}
	// gamma 	:= []float64{18.2, 18.2, 18.2}

	// suksd := []float64{300.0, 0.0, 0.0}
	// phiksd:= []float64{2.6, 0.0, 0.0}
	
	// su11 := []float64{590.0, 590.0, 590.0}
	// phi11:= []float64{35.0, 35.0, 35.0}


	// //////////////////////////////////////////////////////////////////////////
	// /////////////////////////////Birgham//////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////
	// sukmean := []float64{650.0, 650.0, 650.0}
	// phikmean:= []float64{35.0, 35.0, 35.0}
	// gamma 	:= []float64{23.58, 23.58, 23.58}

	// suksd := []float64{300.0, 0.0, 0.0}
	// phiksd:= []float64{2.6, 0.0, 0.0}
	
	// su11 := []float64{650.0, 650.0, 650.0}
	// phi11:= []float64{35.0, 35.0, 35.0}


	// //////////////////////////////////////////////////////////////////////////
	// /////////////////////////////CONE//////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////
	// sukmean := []float64{1444.0, 1444.0, 1444.0}
	// phikmean:= []float64{30.0, 30.0, 30.0}
	// gamma 	:= []float64{20.0, 20.0, 20.0}

	// suksd := []float64{288.0, 0.0, 0.0}
	// phiksd:= []float64{6.0, 0.0, 0.0}
	
	// su11 := []float64{1444.0, 1444.0, 1444.0}
	// phi11:= []float64{30.0, 30.0, 30.0}

	return suksd, phiksd, sukmean, phikmean, gamma, su11, phi11
}
