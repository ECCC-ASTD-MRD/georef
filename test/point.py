# from ctype import Point
import libgeoref.interp_new as libgeoref_new

if __name__ == "__main__":
    point = libgeoref_new.Create(1,2)
    print(point.contents.x)
    # point = Point(10, 20)
    # print(point.x, point.y)