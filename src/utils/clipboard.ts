/**
 * 复制文本到剪贴板：优先 Clipboard API，失败时用隐藏 textarea + execCommand 回退。
 * 两条路径都失败则 reject（调用方可自行处理）。
 */
export async function copyText(text: string): Promise<void> {
  try {
    await navigator.clipboard.writeText(text)
  } catch {
    const textArea = document.createElement('textarea')
    textArea.value = text
    document.body.appendChild(textArea)
    textArea.select()
    document.execCommand('copy')
    document.body.removeChild(textArea)
  }
}
