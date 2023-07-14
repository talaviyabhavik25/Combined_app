import streamlit as st


def main():
    html_temp="""
    <div style="background-color:lightblue;padding:16px">
    <h2 style="color:black";text-align:center> Heater Efficiency Calculator </h2>
    </div>
    """
    st.markdown(html_temp, unsafe_allow_html=True)

    p1 = st.number_input("Enter Oxygen percent %")
    p2 = st.number_input("Enter Stack Temperature in Celsius")
    x= (p1*91.2)/(20.95-2)
    y = (p2*1.8)+32
    efficiency = round((100-((0.0237+(0.000189*x))*(y-80)))*(100/(100+2)),2)

    if st.button("Calculate Efficiency"):
        st.success("Your Heater efficiency is {}".format(efficiency))

if __name__ == '__main__':
    main()



