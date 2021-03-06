////////////////////////////////////////////////////////////////////////
//                 JS-CODE FOR CONTEXTUAL BANDIT                      //
//                       AUTHOR: ERIC SCHULZ                          //
//                    UCL LONDON,  OCTOBER 2014                       //
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//INTIALIZE ARRAYS
////////////////////////////////////////////////////////////////////////
//Number of total trials
var ntrials = 150;
//index for array tracking
var index=0;
//total score
var totalscore=0;
//dummy for trial length
var triallength = new Array(ntrials);
//chosen deck
var chosendeck = new Array(ntrials);
//generated sun tracker
var sunseen = new Array(ntrials);
//generated rain tracker
var rainseen = new Array(ntrials);
//generated temperature tracker
var tempseen = new Array(ntrials);
//received reward tracker
var reward = new Array(ntrials);
//pre-initialize sun, temp, and rain for presentation
var firstgen=generate();
var sun1=firstgen[0];
var temp1 =firstgen[1];
var rain1=firstgen[2];
//k is a dummy to track whether Emeralds or context s currently shown
var k=0;
//chosen deck
var chosen=0;
//producend outcome
var deckproduce=0;


////////////////////////////////////////////////////////////////////////
//CREATE HELPER FUNCTIONS
////////////////////////////////////////////////////////////////////////
//hides page hide and shows page show
function clickStart(hide, show)
{
        document.getElementById(hide).style.display ='none' ;
        document.getElementById(show).style.display ='block';
        window.scrollTo(0,0);        
}

//changes inner HTML of div with ID=x to y
function change (x,y){
    document.getElementById(x).innerHTML=y;
}

//Hides div with id=x
function hide(x){
  document.getElementById(x).style.display='none';
}

//shows div with id=x
function show(x){
  document.getElementById(x).style.display='block';
}

//sets a value at the end to hidden id
function setvalue(x,y){
  document.getElementById(x).value = y;
}

//my norm generation
function myNorm() {
    var x1, x2, rad;
     do {
        x1  = 2 * Math.random() - 1;
        x2  = 2 * Math.random() - 1;
        rad = x1 * x1 + x2 * x2;
    } while(rad >= 1 || rad == 0);
     var c = Math.sqrt(-2 * Math.log(rad) / rad);
     return (x1 * c*5);
};

//Function to randomly shuffle an array and get first element
function generate(){ 
    if (Math.random()>0.5){
    var o=[-1,1,-1];} else{
    var o=[1,-1,1];}
    for(var j, x, i = o.length; i; j = Math.floor(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
    return o;
};

//translate from {1,-1} to {+,-}
function translate(x){
  var y="-";
  if(x===1){
    y="+"
  };
 return(y)
}

hide('submitButton');
////////////////////////////////////////////////////////////////////////
//CREATE EXPERIMENTAL FUNCTIONS
////////////////////////////////////////////////////////////////////////
//deck output generation
function deck(x,y,z,numb) {
 if (numb===1){
//deck1: likes first, dislikes second, indifferent to last, plus noise
    var out=Math.round(50+x*20+y*-20+myNorm())
  } else if(numb===2){
//deck2: dislikes first, indifferent to second, likes last, plus noise
    var out=Math.round(50+x*-20+z*20+myNorm())
  } else if(numb===3){
//deck3: indifferent to first, likes to second, dislikes last, plus noise
    var out=Math.round(50+y*20+z*-20+myNorm())
  } else{
    var out=Math.round(50+myNorm())
  }
  //Show turned on
 k=k+1;
 //return output
 return(out)
};


function showdeck(n){
  //check k
   if (k===0){
   //track chosen deck 
   chosen=n;
   //produce output
  deckproduce=deck(x=sun1,y=temp1,z=rain1, numb=n);
  if (deckproduce>100){deckproduce=100}
  if (deckproduce<0){decckproduce=0}
   //show output on planet
   var whichdeck='deck'+n;
   change(whichdeck,deckproduce);
   //hide context
   ['sun','rain','temp'].forEach(hide);
   //show next button
   show('next');
   //print number of emeralds
   var emeralds='';
  for (i = 0; i < deckproduce; i++) { 
     emeralds=emeralds+'<img src="js/emerald.JPG" width="50" height="50">'
    }
  change('emeralds', emeralds);
  //feedback
  var present='You have mined '+deckproduce+' Emeralds';
  change("outcome",present);
 return(chosen,deckproduce);
 }
}

//overal trial function
function mynexttrial(){
 if (triallength.length > 0) {
 //track trial number
  triallength.shift();
  //track produced output
  sunseen[index]=sun1;rainseen[index]=rain1;tempseen[index]=temp1;
  chosendeck[index]=chosen;
  reward[index]=deckproduce;
  //generate new output
 var gen=generate();
 sun1=gen[0];
 temp1 =gen[1];
 rain1=gen[2];
  //reset show dummy
  k=0;
  //reset button description
  change('deck1','PLANET 1');change('deck2','PLANET 2');
  change('deck3','PLANET 3');change('deck4','PLANET 4');
  //show new context
  var insertsun ='MERCURY: '  +translate(sun1);
  var insertrain ='KRYPTON: ' +translate(rain1);
  var inserttemp ='NOBELIUM: '+translate(temp1);
  //insert new context
  change('sun',insertsun);change('rain', insertrain);change('temp',inserttemp);
  ['sun','rain','temp'].forEach(show);
  //hide next button
  hide('next');
  //keep track for index used to assign array content
  index=index+1;
  //show remaining number of trials
  var insert ='Number of trials left: '+(ntrials-index);
  change("remain",insert);
  //show total score
  totalscore =totalscore+deckproduce;
  var insertscore ='Number of Emeralds harvested: '+totalscore;
  change('score',insertscore);
  //Last trial:
  if((ntrials-index)===0){
   //show "Go to next page button"
   change('nexttrial','Go to next page');
   //Show info that trials are done
   change('finaltext','You have used up all of your trials. Please click on "Go to next page" to continue.');
  }
 }
 else {
   //go to next page
   clickStart('page4','page5');
 }
}

//Final submit function
function mysubmit(){
  //claculate number of mined emeralds overall
 var presenttotal='You have mined a total amount of '+totalscore+' Emeralds.';
 //calculate money earned
 var money =Math.round(100*(0.5+totalscore/(150*100)*0.5))/100;
 var presentmoney='This equals a total reward of $ '+money+'.';
 //show score and money
 change('result',presenttotal); change('money',presentmoney);
 //save all created values
 setvalue('sunseen',sunseen); setvalue('rainseen',rainseen); setvalue('tempseen',tempseen); 
 setvalue('chosendeck',chosendeck); setvalue('reward', reward); setvalue('total', totalscore); setvalue('money', money);
 //Go to final "Thank you"-page
 clickStart('page5','page6');
}
////////////////////////////////////////////////////////////////////////
