package Part

import (
	// "fmt"
	"io/ioutil"
	"strconv"
	"os/exec"
	"bufio"
	"log"
	"strings"
	"os"
)

func Partmain(x, y, z, R, alpha float64) (score float64, err bool){
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////SINGLE SOIL MATERIAL/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   41.65   15   18.82"

	// FileName := "StHelens"
	// ASCFile := "sthel_res100mDEM.asc"
	// material := "1   1000.0   40.0   24.0"
	
	// FileName := "Cone"
	// ASCFile := "cone1000DEM.asc"
	// material := "1   17630   40   21"
	// material := "1   1444   30   20"
	// material := "1   172.9   20   19"
	
	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\noutput directory\nOutput/")

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////MULTIPLE MATERIAL/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material1 := "1   29.4   12   18.82"
	// material2 := "2   9.8   5   18.82"
	// material3 := "3   294   40   18.82" 
	// LayerFile := "emb20layer"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nnmat\n3\nlnum   cee   phi   gamt\n"+material1+"\n"+material2+"\n"+material3+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   5.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\nlayer file\n"+FileName+"/input/"+LayerFile+"\noutput directory\nOutput/")

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////3D MULTIPLE MATERIAL DATA////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   -1   -1   -1"
	// MaterialFile := "emb20mat3D.txt"

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
	FileName := "Cone"
	ASCFile := "cone1000DEM.asc"
	material := "1   1444   30   20   0.2"

	d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nru\nnmat\n1\nlnum   cee   phi   gamt   ru\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\noutput directory\nOutput/")

	err1 := ioutil.WriteFile("test.scp", d1, 0644)
	Check(err1)

	//exe.go --> 
	exec.Command("/home/yewintun/mygo/src/scoops3d/fortrancode/scoops3d.exe").Run()
	
	//read.go --> read test_out.txt -->FOS
	file, err3 := os.Open("/home/yewintun/mygo/src/scoops3d/Output/test_out.txt")

	if err3 != nil {
		log.Fatal(err3)
	}

	defer file.Close()

	important, important2 := "1st", "2nd"
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
		}
	}
	
	important3 := strings.Fields(important2)
	// fmt.Println(important3)
	if important3[0] == "2nd" {
		score = 1000.0
	}else{
		score, _ =(strconv.ParseFloat(important3[5], 64))
		if score == 0.0{
			score = 1000.0
		}
	}

	if err3 := scanner.Err(); err3 != nil {
		log.Fatal(err3)
	}

	if err3 := scanner.Err(); err3 != nil {
		log.Fatal(err3)
	}
	return score, err
}

func Check(e error) {
	if e != nil {
		panic(e)
	}
}

