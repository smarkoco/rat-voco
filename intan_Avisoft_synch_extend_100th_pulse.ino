
// https://www.instructables.com/Arduino-Timer-Interrupts/

const byte inputPin = 2; // start recording pin
const byte intanPin = 11; // intan output
const byte outputPin = 12; // camera output

const byte modePin = A0; // mode select
const byte ledPin = LED_BUILTIN; // led pin

int mode = 0; // could be used to set different frequencies

#define MODE_0 0
#define MODE_1 1

volatile boolean record_on = 0;
// volatile boolean is_first_call = 1;
volatile boolean toggle1 = 0;
volatile uint8_t counter_for_100th_pulse = 0;

void setup() {
  // put your setup code here, to run once:

  pinMode(inputPin, INPUT);
  pinMode(modePin, INPUT_PULLUP);
  pinMode(outputPin, OUTPUT);
  pinMode(ledPin, OUTPUT);
  
  cli();

  //set timer1 interrupt at desired freq
  TCCR1A = 0; // set entire TCCR1A register to 0
  TCCR1B = 0; // same for TCCR1B
  TCNT1  = 0; //initialize counter value to 0
  // set compare match register for 100hz increments (2*freq)
  OCR1A = 9999; // = (16*10^6) / (2*freq*8) - 1 (must be <65536)
  // turn on CTC mode
  TCCR1B |= (1 << WGM12);
  // Set CS11 bits for 8 prescaler
  TCCR1B |= (1 << CS11);  
  // enable timer compare interrupt
  TIMSK1 |= (1 << OCIE1A);

  attachInterrupt(digitalPinToInterrupt(inputPin), start_timer, RISING);

  sei(); //allow interrupts

  //Serial.begin(38400);
}

ISR(TIMER1_COMPA_vect){ 
  //generates pulse wave of timer freq/2 (takes two cycles for full wave- toggle high then toggle low)
  if (record_on){
    if (counter_for_100th_pulse<=197){
      if (toggle1){
          digitalWrite(outputPin,HIGH);
          toggle1 = 0;
      }
      else
      {
          digitalWrite(outputPin,LOW);
          toggle1 = 1;
      }
      counter_for_100th_pulse++;
    }
    // for 99th pulse, keep it HIGH, i.e., do not change digitalWrite() value; also increment counter
    if (counter_for_100th_pulse>197 && counter_for_100th_pulse<201){
      counter_for_100th_pulse++;
    }
    // at next half pulse increment, keep LOW, and reset counter, then next pulse will be HIGH
    if (counter_for_100th_pulse>=201){
      counter_for_100th_pulse = 0;
    }
  }
  else{
    digitalWrite(outputPin,LOW);
  }
}

void start_timer()
{
    // reset timer count
    TCNT1 = 0; 
    // set timer and intan output to HIGH
    digitalWrite(outputPin,HIGH);
    digitalWrite(intanPin,HIGH);

    // timer output starts at HIGH, so next value must be LOW (0) 
    toggle1 = 0;
    
    // reset counter value at new recording start
    counter_for_100th_pulse = 0;
}


void loop() {

  mode = 1-digitalRead(modePin);
  digitalWrite(ledPin, mode);

  record_on = digitalRead(inputPin);
  digitalWrite(intanPin, record_on);

  // ensure counter is reset after the 99th pulse
//  if (counter_for_100th_pulse > 200){
//    counter_for_100th_pulse = 0;
//    }
}
