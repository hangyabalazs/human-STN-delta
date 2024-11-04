/* Timestamping some input button presses
   Andras Szell
*/

// debug: print "AAA" as long as button is pressed
// 1 to enable, 0 when no need for debug - corresponding lines are to be removed from production code :)
#define DEBUGKEYS 0

const int pinCnt = 5;
const long debounceTimeMs = 10; // the debounce time, increase if the output flickers
// but then key presses shorter than this will be filtered out
const int inPins[] = {10, 11, 12, 13, 9};    // number of input pins listed in this array should match pinCnt

long debounceTimeoutMs[pinCnt]; // 0 if no debouncing on the given pin, future timeout value if debouncing is in progress
bool prevState[pinCnt];         // previous button state
int  pressCnt[pinCnt];          // number of button presses
bool prevWhite;         // previous light state

// define photosensor keys
int photoRPin = 0;     //Set the photodiode pin
int triggerPin = 19;
int minLight;          //Used to calibrate the readings
int maxLight;          //Used to calibrate the readings
int lightLevel;
int prevLight;
int UpLight;
int LowLight;

void setup()
{
  Serial.begin(115200); // 115200 baud -> better timing of serial messages if necessary
  Serial.println("Debouncer started");

  // safety check :)
  if (sizeof(inPins) / sizeof(int) != pinCnt) {
    Serial.println("Error: wrong number of input pins in pinCnt!");
  }

  // set all pins listed in inPins to input
  for (int i = 0; i < pinCnt; i++) {
    pinMode(inPins[i], INPUT);
  }

  // set up the starting light level limits
  lightLevel = analogRead(photoRPin);
  minLight = 300;
  maxLight = lightLevel;
  prevLight = maxLight;
  prevWhite = HIGH;
  pinMode(inPins[4], OUTPUT);
  pinMode(triggerPin, OUTPUT);
  digitalWrite(triggerPin, LOW);
}

int val = 0;

void loop()
{
  long now = millis();

  // Light measurements
  lightLevel = analogRead(photoRPin);
  if (minLight > lightLevel) {
    minLight = lightLevel;
  }
  if (maxLight < lightLevel) {
    maxLight = lightLevel;
  }
  UpLight = maxLight - 200;
  LowLight = minLight + 200;

  // Set digital pin based on the analog input from the sensor
  if (lightLevel < UpLight) {   // Dark rectangle in the top left corner (Cue)
    //if (prevWhite == HIGH) {
      digitalWrite(inPins[4], HIGH);
      digitalWrite(triggerPin, HIGH);
      //Serial.println(lightLevel);
      //Serial.println("fekete");
    //}
  } else if (lightLevel > LowLight) {
    digitalWrite(inPins[4], LOW);   // White background
    digitalWrite(triggerPin, LOW);
    //prevWhite = HIGH;
    //Serial.println(lightLevel);
    //Serial.println("feher");
  }


  for (int i = 0; i < pinCnt; i++) {
#if DEBUGKEYS
    delay(1);
#endif

    // check pins one by one, but only if not in debounce timeout's 50 ms window

    if (debounceTimeoutMs[i] == 0) {
      // nothing was pressed recently - wait for button press and always update previous state

      if (digitalRead(inPins[i]) == HIGH) {
        if (prevState[i] == LOW) {  // rising edge: button is really pressed -> setup timeout
          debounceTimeoutMs[i] = now + debounceTimeMs;
        }
        prevState[i] = HIGH;
#if DEBUGKEYS
        Serial.print((char)('A' + i));
#endif
      } else {
        prevState[i] = LOW;
      }

    } else if (debounceTimeoutMs[i] < now) {

      // debounce time is over
      if (digitalRead(inPins[i]) == HIGH) {
        //  we assume a real button press if it is still HIGH
        pressCnt[i]++;
#if DEBUGKEYS
        Serial.print(' ');
#endif
        Serial.print(i); Serial.print(" at ");
        Serial.print(debounceTimeoutMs[i] - debounceTimeMs);
        Serial.print(" press ");
        Serial.print(pressCnt[i]);
        Serial.println("");
        //prevWhite = LOW;
        //digitalWrite(inPins[4], LOW);
        //digitalWrite(triggerPin, LOW);
        //Serial.println(lightLevel);
        //Serial.println("meg mindig fekete");
      } else {
#if DEBUGKEYS
        Serial.print("bounce "); Serial.println(i);
#endif
        prevState[i] = LOW;
      }
      debounceTimeoutMs[i] = 0;
    }
    // else: debounceTimeout is not over, keep waiting...
  }
}


