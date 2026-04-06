from app import SpectroApp
import tkinter as tk
from PIL import Image, ImageTk
import os

def show_splash_screen():
    """Display splash screen with SciGlob logo with fade-in effect."""
    splash = tk.Tk()
    splash.title("SciGlob")
    splash.overrideredirect(True)  # Remove window decorations
    
    # Get screen dimensions
    screen_width = splash.winfo_screenwidth()
    screen_height = splash.winfo_screenheight()
    
    # Set splash window size
    splash_width = 700
    splash_height = 450
    
    # Calculate position to center the splash screen
    x = (screen_width - splash_width) // 2
    y = (screen_height - splash_height) // 2
    
    splash.geometry(f"{splash_width}x{splash_height}+{x}+{y}")
    
    # Set dark background to match professional look
    bg_color = '#1a1a1a'  # Dark gray background
    splash.configure(bg=bg_color)
    
    # Start with fully transparent
    splash.attributes('-alpha', 0.0)
    
    # Load and display logo
    try:
        logo_path = os.path.join(os.path.dirname(__file__), "sciglob_logoRGB.png")
        img = Image.open(logo_path)
        
        # Resize image to fit splash screen while maintaining aspect ratio
        img.thumbnail((600, 380), Image.Resampling.LANCZOS)
        photo = ImageTk.PhotoImage(img)
        
        # Create label with no borders or outlines
        logo_label = tk.Label(splash, image=photo, 
                             bg=bg_color,
                             borderwidth=0, 
                             highlightthickness=0)
        logo_label.image = photo  # Keep a reference
        logo_label.pack(expand=True)
        
    except Exception as e:
        # Fallback if logo can't be loaded
        label = tk.Label(splash, text="SciGlob\nSpectrometer Characterization System", 
                        font=("Arial", 20, "bold"), bg=bg_color, fg='#00a8e8',
                        borderwidth=0, highlightthickness=0)
        label.pack(expand=True)
        print(f"Could not load logo: {e}")
    
    # Update window
    splash.update()
    
    # Fade-in animation
    def fade_in(alpha=0.0):
        if alpha < 1.0:
            alpha += 0.05  # Increment alpha
            splash.attributes('-alpha', alpha)
            splash.after(30, lambda: fade_in(alpha))  # Call again after 30ms
    
    # Start fade-in effect
    fade_in()
    
    # Auto-close after 2.5 seconds
    splash.after(2500, splash.destroy)
    splash.mainloop()

if __name__ == '__main__':
    # Show splash screen
    show_splash_screen()
    
    # Create and run main application
    app = SpectroApp()
    app.mainloop()
