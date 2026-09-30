const std = @import("std");
const Io = std.Io;

const zmath = @import("zmath");

pub fn main(init: std.process.Init) !void {
    // This is appropriate for anything that lives as long as the process.
    const arena: std.mem.Allocator = init.arena.allocator();

    // Accessing command line arguments:
    const args = try init.minimal.args.toSlice(arena);
    for (args) |arg| {
        std.log.info("arg: {s}", .{arg});
    }

    // In order to do I/O operations need an `Io` instance.
    const io = init.io;

    // Stdout is for the actual output of your application, for example if you
    // are implementing gzip, then only the compressed bytes should be sent to
    // stdout, not any debugging messages.
    var stdout_buffer: [1024]u8 = undefined;
    var stdout_file_writer: Io.File.Writer = .init(.stdout(), io, &stdout_buffer);
    const stdout_writer = &stdout_file_writer.interface;

    // Report the result bits, not a rounded decimal, because bit identity to the
    // scalar kernel is the contract this library sells.
    const result: f64 = zmath.exp2m1(@as(f64, 1.0));
    try stdout_writer.print("exp2m1(1.0) = 0x{x}\n", .{@as(u64, @bitCast(result))});

    try stdout_writer.flush(); // Don't forget to flush!
}
