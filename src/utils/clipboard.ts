/**
 * 复制文本到剪贴板：优先 Clipboard API，失败时用隐藏 textarea + execCommand 回退。
 * 成功返回 true，两条路径都失败返回 false。
 */
export async function copyText(text: string): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text)
    return true
  } catch {
    const activeElement = document.activeElement
    const selection = document.getSelection()
    const savedRange =
      selection && selection.rangeCount > 0 ? selection.getRangeAt(0).cloneRange() : null

    const textArea = document.createElement('textarea')
    textArea.value = text
    textArea.setAttribute('readonly', '')
    textArea.style.position = 'fixed'
    textArea.style.top = '0'
    textArea.style.left = '0'
    textArea.style.opacity = '0'
    document.body.appendChild(textArea)
    textArea.focus()
    textArea.select()
    textArea.setSelectionRange(0, textArea.value.length)

    let copied = false
    try {
      copied = document.execCommand('copy')
    } catch {
      copied = false
    }

    document.body.removeChild(textArea)

    if (activeElement instanceof HTMLElement) activeElement.focus()
    if (savedRange && savedRange.startContainer.isConnected) {
      const restored = document.getSelection()
      restored?.removeAllRanges()
      restored?.addRange(savedRange)
    }
    return copied
  }
}
