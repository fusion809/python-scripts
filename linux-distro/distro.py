time = input("How long are you willing to spend setting up your Linux system? This time should be in hours, if you want help type 'help' ")
if time=='help':
    print("Your Linux system will not be able to do everything you want it to immediately, you will need to configure it according to your own preferences. So this question is essentially how much time you, personally, are willing to spend configuring your system and making it suit your needs. If you want a system that should require minimal configuration then give an answer that is less than 1. If you want a system you have to build from the ground-up (potentially taking > 5 hours), you should choose a larger number.")
    time = input("How long are you willing to spend setting up your Linux system? This time should be in hours, if you want help type 'help' ")

###
patience = input("How patient are you, from one to ten? With ten being you are a saint, with respect to this particular virtue, and one being you have little patience. ")

foss=input("Where do you land on the FOSS spectrum, from 1 to 5? 1 being FOSS is not important, proprietary software is fine! 5 being it is either FOSS or nothing! Enter 'help' for further info. ")

if foss=='help':
    print("FOSS or free and open-source software is software that is licensed such that its source code and pre-built forms can be freely and openly shared, modified, redistributed, etc. without legal restriction.")
