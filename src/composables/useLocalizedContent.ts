import { ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'

/**
 * 内容页统一的数据加载模式：
 * 首次加载 + locale 变更时重载，异常回退到 defaultValue。
 * 覆盖 Home/Category/Resource/About/PostList/NavigationTree 的重复样板。
 */
export function useLocalizedContent<T>(load: () => T, defaultValue: T) {
  const data = ref<T>(defaultValue)
  const { locale } = useI18n()

  const reload = () => {
    try {
      data.value = load()
    } catch (err) {
      console.error('Failed to load localized content:', err)
      data.value = defaultValue
    }
  }

  reload()
  watch(locale, reload)

  return { data, reload }
}
